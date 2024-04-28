package progistar.scan.run;

import java.io.File;
import java.util.ArrayList;
import java.util.List;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import progistar.scan.data.BAMSRecord;
import progistar.scan.data.Constants;
import progistar.scan.data.LocTable;

public class Task {

	public int taskIdx = -1;
	public int type = Constants.TYPE_MAPPED_TASK;
	public ArrayList<BAMSRecord> records = new ArrayList<BAMSRecord>();
	
	// only available for FullMode.
	public LocTable locTable = new LocTable();
	public String chrName;
	public int start;
	public int end;
	
	public String getTaskInfo () {
		String typeStr = "*mapped";
		if(type == Constants.TYPE_UNMAPPED_TASK) {
			typeStr = "*unmapped";
		}else if(type == Constants.TYPE_DISCOVERY_TASK) {
			typeStr = "discovery";
		}
		return "Task"+taskIdx+" has "+records.size()+" "+typeStr+" records";
	}
	
	public Task(int type) {
		this.type = type;
	}
	
	
	private static ArrayList<Task> getChromosomeLevelTasks (ArrayList<BAMSRecord> records, String chrName, int start, int end) {
		ArrayList<Task> tasks = new ArrayList<Task>();
		
		int size = end - start + 1;
		int interval = size / Scan.threadNum + 1;
		
		int startInterval = 1;
		int endInterval = 1;
		
		boolean isEndOfLength = false;
		for(int i=0; i<Scan.threadNum; i++) {
			endInterval = startInterval + interval - 1;
			
			if(endInterval >= end) {
				endInterval = end;
				isEndOfLength = true;
			}
			
			Task task = new Task(Constants.TYPE_DISCOVERY_TASK);
			task.records = records;
			task.chrName = chrName;
			task.start = startInterval;
			task.end = endInterval;
			tasks.add(task);
			
			startInterval = endInterval + 1;
			
			if(isEndOfLength) {
				break;
			}
		}
		
		return tasks;
	}
	
	public static ArrayList<Task> getFullTask (ArrayList<BAMSRecord> records) {
		ArrayList<Task> tasks = new ArrayList<Task>();
		
		File file = new File(Scan.bamFile.getAbsolutePath());
		
		try (SamReader samReader = SamReaderFactory.makeDefault().open(file)) {
			System.out.println(samReader.getFileHeader().getSequenceDictionary().getSequences().get(0).getSequenceLength());
			List<SAMSequenceRecord> chromosomes = samReader.getFileHeader().getSequenceDictionary().getSequences();
			for(SAMSequenceRecord chromosome : chromosomes) {
				System.out.println(chromosome.getSAMString());
				String chrName = chromosome.getSequenceName();
				int start = chromosome.getStart();
				int end = chromosome.getEnd();
				
				tasks.addAll(getChromosomeLevelTasks(records, chrName, start, end));
			}
			
			// for unmapped reads
			SAMRecordIterator unmappedIter = samReader.queryUnmapped();
			int size = 0;
			while(unmappedIter.hasNext()) {
				SAMRecord samRecord = unmappedIter.next();
				if(Scan.unmmapedMarker == null) {
					Scan.unmmapedMarker = samRecord.getReferenceName();
				}
				size ++;
			}
			if(Scan.unmmapedMarker != null) {
				System.out.println("@SQ\t"+Scan.unmmapedMarker+"\tLN:"+size);
				tasks.addAll(getChromosomeLevelTasks(records, Scan.unmmapedMarker, 1, size));
			}
			
			// assign idx
			for(int i=0; i<tasks.size(); i++) {
				tasks.get(i).taskIdx = (i+1);
			}
			
			unmappedIter.close();
		} catch(Exception e) {
			e.printStackTrace();
			System.exit(1);
		}
		
		
		return tasks;
	}
	
	public static ArrayList<Task> divideTasks (ArrayList<BAMSRecord> records, int chunkSize) {
		System.out.println("Prepare tasks with chunk size = "+chunkSize);
		
		ArrayList<Task> tasks = new ArrayList<Task>();
		
		ArrayList<BAMSRecord> mappedRecords = new ArrayList<BAMSRecord>();
		ArrayList<BAMSRecord> unmappedRecords = new ArrayList<BAMSRecord>();
		records.forEach(BAMSRecord -> {
			if(BAMSRecord.chr.equalsIgnoreCase("*")) {
				// unmapped classes
				unmappedRecords.add(BAMSRecord);
			} else {
				// mapped classes
				mappedRecords.add(BAMSRecord);
			}
		});
		
		int mappedSize = mappedRecords.size();
		int unmappedSize = unmappedRecords.size();
		
		System.out.println("Mapped records: "+mappedSize);
		System.out.println("=> "+(mappedSize/chunkSize+1) +" tasks");
		System.out.println("Unmapped records: "+unmappedSize);
		System.out.println("=> "+(unmappedSize/chunkSize+1) +" tasks");
		
		// generate unmapped tasks
		int sIdx = 0;
		while(sIdx < unmappedSize) {
			int eIdx = sIdx + chunkSize > unmappedSize ? unmappedSize : sIdx + chunkSize;
			Task task = new Task(Constants.TYPE_UNMAPPED_TASK);
			for(int i=sIdx; i<eIdx; i++) {
				task.records.add(unmappedRecords.get(i));
			}
			tasks.add(task);
			task.taskIdx = tasks.size();
			sIdx = eIdx;
		}
		
		// generate mapped tasks
		sIdx = 0;
		while(sIdx < mappedSize) {
			int eIdx = sIdx + chunkSize > mappedSize ? mappedSize : sIdx + chunkSize;
			Task task = new Task(Constants.TYPE_MAPPED_TASK);
			for(int i=sIdx; i<eIdx; i++) {
				task.records.add(mappedRecords.get(i));
			}
			tasks.add(task);
			task.taskIdx = tasks.size();
			sIdx = eIdx;
		}
		
		System.out.println("Task list");
		for(Task task : tasks) {
			System.out.println(task.getTaskInfo());
		}
		
		return tasks;
	}
}
