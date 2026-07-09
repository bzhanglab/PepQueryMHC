package progistar.scan.function;

import java.io.File;
import java.util.ArrayList;
import java.util.LinkedList;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

import htsjdk.samtools.fastq.FastqReader;
import htsjdk.samtools.fastq.FastqRecord;
import progistar.scan.data.BarcodeTable;
import progistar.scan.data.Constants;
import progistar.scan.data.LocTable;
import progistar.scan.data.Parameters;
import progistar.scan.run.Task;

public class FASTQModeRun extends Mode {

	public static void runFASTQMode (Task task) {
		if(task.type == Constants.TYPE_FASTQ_MODE_TASK) {
			if(Parameters.verbose) {
				System.out.println(task.chrName+":"+task.start+"-"+task.end);
			}
			scanReads(task);
		}
	}
	
	private static void scanReads (Task task) {
		long startTime = System.currentTimeMillis();
		File file = null;
		int threadNum = Parameters.threadNum+1;
		
		// single-end
		if(task.start == 0) {
			file = new File(Parameters.fastq0File.getAbsolutePath());
		} 
		// first fastq
		else if(task.start == 1) {
			file = new File(Parameters.fastq1File.getAbsolutePath());
			threadNum /= 2;
		} 
		// second fastq
		else if(task.start == 2) {
			file = new File(Parameters.fastq2File.getAbsolutePath());
			threadNum /= 2;
		}
		
		
		ExecutorService executorService = Executors.newFixedThreadPool(threadNum);
		// Pending chunk futures. Submission and result-merging are decoupled so that up to
		// "threadNum" chunks can be processed concurrently instead of blocking on each one.
		// Only the reading thread ever calls task.locTable.putLocation(...), so no
		// synchronization is needed around the (non-thread-safe) HashMap-backed LocTable.
		LinkedList<Future<LocTable>> pendingChunks = new LinkedList<Future<LocTable>>();
		try(FastqReader reader = new FastqReader(file)) {
			int maxSize = 10000;
			int curSize = 0;

			ArrayList<FastqRecord> records = new ArrayList<FastqRecord>();
			for(FastqRecord fastqRecord : reader) {
				String barcodeId = BarcodeTable.getBarcodeFromFASTQ(fastqRecord);

	    		// increase processed reads
	        	Double pReads = task.processedReads.get(barcodeId);
	        	if(pReads == null) {
	        		pReads = .0;
	        	}
	        	pReads++;
	        	task.processedReads.put(barcodeId, pReads);

	        	//
				records.add(fastqRecord);
				curSize++;

				if(curSize == maxSize || !reader.hasNext()) {
					final ArrayList<FastqRecord> copies = records;

					// note that if you process something inside "find" function,
					// you are care about concurrence, conflicts.
					pendingChunks.add(executorService.submit(() -> find(copies, Task.allTrie, task)));
					records = new ArrayList<FastqRecord>();
					curSize = 0;

					// keep at most "threadNum" chunks in flight so more than one can run
					// concurrently, while still bounding memory usage.
					while(pendingChunks.size() >= threadNum) {
						mergeChunkResult(task, pendingChunks.poll().get());
					}
				}
			}

			// drain remaining in-flight chunks
			while(!pendingChunks.isEmpty()) {
				mergeChunkResult(task, pendingChunks.poll().get());
			}
		} catch (Exception e) {
			e.printStackTrace();
		} finally {
			executorService.shutdown();
		}

		while(!executorService.isTerminated()) {
			try {
				Thread.sleep(1000);
			} catch (InterruptedException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
		long endTime = System.currentTimeMillis();
		if(Parameters.verbose) {
			System.out.println("Task"+task.taskIdx+" "+(endTime-startTime)/1000+" sec");
		}
	}

	private static void mergeChunkResult (Task task, LocTable locTable) {
		locTable.table.forEach((key, subTable)->{
			subTable.forEach((key_, linfo)->{
				task.locTable.putLocation(linfo);
			});
		});
	}
}
