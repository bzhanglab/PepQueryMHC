package progistar.scan.function;

import java.io.File;
import java.util.ArrayList;
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
		int threadNum = Parameters.threadNum;
		
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
		try(FastqReader reader = new FastqReader(file)) {
			int maxSize = 10000;
			int curSize = 0;
			int debugTotal = 1000000;
			
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
					Future<LocTable> future = executorService.submit(() -> find(copies, Task.allTrie, task));
					records = new ArrayList<FastqRecord>();
					curSize = 0;
					
					LocTable locTable = future.get();
					locTable.table.forEach((key, subTable)->{
						subTable.forEach((key_, linfo)->{
							task.locTable.putLocation(linfo);
						});
					});
				}
				
				debugTotal --;
				if(debugTotal == 0) {
					break;
				}
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
}
