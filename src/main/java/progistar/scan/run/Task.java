package progistar.scan.run;

import java.util.ArrayList;
import java.util.Hashtable;
import java.util.LinkedList;

import progistar.scan.data.BAMSRecord;
import progistar.scan.data.Constants;

public class Task {

	public int taskIdx = -1;
	public int type = Constants.TYPE_MAPPED_TASK;
	public ArrayList<BAMSRecord> records;
	
	public String getTaskInfo () {
		return "Task"+taskIdx+" has "+records.size()+" "+type+" records";
	}
	
	public Task(int type) {
		this.type = type;
	}
	
	public static LinkedList<Task> divideTasks (ArrayList<BAMSRecord> records, int chunkSize) {
		System.out.println("Prepare tasks with chunk size = "+chunkSize);
		
		LinkedList<Task> tasks = new LinkedList<Task>();
		Hashtable<String, BAMSRecord> uniqueRecords = new Hashtable<String, BAMSRecord>();
		
		records.forEach((record)->{
			String key = record.getKey();
			uniqueRecords.put(key, record);
		});
		
		System.out.println("Unique records: "+uniqueRecords.size());
		
		ArrayList<BAMSRecord> mappedRecords = new ArrayList<BAMSRecord>();
		ArrayList<BAMSRecord> unmappedRecords = new ArrayList<BAMSRecord>();
		uniqueRecords.forEach((key, BAMSRecord)->{
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
			task.records = (ArrayList<BAMSRecord>) unmappedRecords.subList(sIdx, eIdx);
			tasks.add(task);
			task.taskIdx = tasks.size();
			sIdx = eIdx;
		}
		
		// generate mapped tasks
		sIdx = 0;
		while(sIdx < mappedSize) {
			int eIdx = sIdx + chunkSize > mappedSize ? mappedSize : sIdx + chunkSize;
			Task task = new Task(Constants.TYPE_MAPPED_TASK);
			task.records = (ArrayList<BAMSRecord>) mappedRecords.subList(sIdx, eIdx);
			tasks.add(task);
			task.taskIdx = tasks.size();
			sIdx = eIdx;
		}
		
		System.out.println("Task list");
		for(Task task : tasks) {
			String type = "mapped";
			if(task.type == Constants.TYPE_UNMAPPED_TASK) {
				type = "unmapped";
			}
			System.out.println(task.getTaskInfo());
		}
		
		return tasks;
	}
}
