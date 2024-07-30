package progistar.scan.run;

import java.io.File;
import java.util.ArrayList;
import java.util.List;

import org.ahocorasick.trie.Trie;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import progistar.scan.data.Constants;
import progistar.scan.data.LocTable;
import progistar.scan.data.SequenceRecord;

public class Task implements Comparable<Task> {

	public int taskIdx = -1;
	public int type = Constants.TYPE_TARGET_MODE_MAPPED_TASK;
	public int readType = Constants.MAPPED_READS;
	public double processedReads = 0;
	public ArrayList<SequenceRecord> records = new ArrayList<SequenceRecord>();
	
	// only available for FullMode.
	public static Trie allTrie = null;
	public LocTable locTable = new LocTable();
	public String chrName;
	public int start;
	public int end;
	
	// estimated memory usage
	public long peakMemory = 0;
	
	
	public String getTaskInfo () {
		String typeStr = "target mode (mapped)";
		if(type == Constants.TYPE_TARGET_MODE_UNMAPPED_TASK) {
			typeStr = "target mode (unmapped)";
		} else if(type == Constants.TYPE_SCAN_MODE_TASK) {
			typeStr = "scan mode";
		} else if(type == Constants.TYPE_TARGET_MODE_LIBRARY_ESTIMATION_TASK) {
			typeStr = "target mode (libary estimation)";
		}
		return typeStr+": Task"+taskIdx+" has "+records.size();
	}
	
	public Task(int type) {
		this.type = type;
	}
	
	
	private static ArrayList<Task> getChromosomeLevelTasks (ArrayList<SequenceRecord> records, String chrName, int start, int end) {
		ArrayList<Task> tasks = new ArrayList<Task>();
		
		int divider = Scan.threadNum * Scan.chunkSize;
		// if unmapped reads?
		// only number of threads will be partitioned
		if(chrName.equalsIgnoreCase(Constants.NULL)) {
			divider = Scan.threadNum;
		}
		int size = end - start + 1;
		int interval = Math.max( (size / divider) + 1, 100000);
		
		int startInterval = 1;
		int endInterval = 1;
		
		boolean isEndOfLength = false;
		while(true) {
			endInterval = startInterval + interval - 1;
			
			if(endInterval >= end) {
				endInterval = end;
				isEndOfLength = true;
			}
			
			Task task = new Task(Constants.TYPE_SCAN_MODE_TASK);
			// if the chrName equals to Constans.NULL => it is belonged to an unmapped read group
			if(chrName.equalsIgnoreCase(Constants.NULL)) {
				task.readType = Constants.UNMAPPED_READS;
			} else {
				task.readType = Constants.MAPPED_READS;
			}
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
	
	public static ArrayList<Task> getScanModeTasks (ArrayList<SequenceRecord> records) {
		ArrayList<Task> tasks = new ArrayList<Task>();
		
		File file = new File(Scan.bamFile.getAbsolutePath());
		// build global trie
		System.out.println("Build Trie");
		Task.allTrie = SequenceRecord.getTrie(records);
		System.out.println("Complete building Trie");
		try (SamReader samReader = SamReaderFactory.makeDefault().open(file)) {
			//System.out.println(samReader.getFileHeader().getSequenceDictionary().getSequences().get(0).getSequenceLength());
			List<SAMSequenceRecord> chromosomes = samReader.getFileHeader().getSequenceDictionary().getSequences();
			for(SAMSequenceRecord chromosome : chromosomes) {
				// System.out.println(chromosome.getSAMString());
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
				// System.out.println("@SQ\t"+Scan.unmmapedMarker+"\tLN:"+size);
				tasks.addAll(getChromosomeLevelTasks(records, Constants.NULL, 1, size));
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
	
	public static ArrayList<Task> getTargetModeTasks (ArrayList<SequenceRecord> records, int chunkSize) {
		System.out.println("Prepare tasks with chunk size = "+chunkSize);
		
		ArrayList<Task> tasks = new ArrayList<Task>();
		
		ArrayList<SequenceRecord> mappedRecords = new ArrayList<SequenceRecord>();
		ArrayList<SequenceRecord> unmappedRecords = new ArrayList<SequenceRecord>();
		records.forEach(BAMSRecord -> {
			if(BAMSRecord.chr.equalsIgnoreCase(Constants.NULL)) {
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
		//System.out.println("=> "+(unmappedSize/chunkSize+1) +" tasks");
		System.out.println("=> "+(1) +" task");
		
		// generate unmapped tasks
		int sIdx = 0;
		while(sIdx < unmappedSize) {
			/*
			int eIdx = sIdx + chunkSize > unmappedSize ? unmappedSize : sIdx + chunkSize;
			Task task = new Task(Constants.TYPE_UNMAPPED_TASK);
			for(int i=sIdx; i<eIdx; i++) {
				task.records.add(unmappedRecords.get(i));
			}
			tasks.add(task);
			task.taskIdx = tasks.size();
			sIdx = eIdx;*/
			Task task = new Task(Constants.TYPE_TARGET_MODE_UNMAPPED_TASK);
			task.readType = Constants.UNMAPPED_READS;
			task.records = unmappedRecords;
			tasks.add(task);
			task.taskIdx = tasks.size();
			sIdx = unmappedSize;
		}
		
		// generate mapped tasks
		sIdx = 0;
		while(sIdx < mappedSize) {
			int eIdx = sIdx + chunkSize > mappedSize ? mappedSize : sIdx + chunkSize;
			Task task = new Task(Constants.TYPE_TARGET_MODE_MAPPED_TASK);
			task.readType = Constants.MAPPED_READS;
			
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
	
	public static ArrayList<Task> getLibSizeTask (ArrayList<SequenceRecord> records) {
		ArrayList<Task> tasks = new ArrayList<Task>();
		
		File file = new File(Scan.bamFile.getAbsolutePath());
		Task.allTrie = null;
		
		try (SamReader samReader = SamReaderFactory.makeDefault().open(file)) {
			System.out.println(samReader.getFileHeader().getSequenceDictionary().getSequences().get(0).getSequenceLength());
			List<SAMSequenceRecord> chromosomes = samReader.getFileHeader().getSequenceDictionary().getSequences();
			for(SAMSequenceRecord chromosome : chromosomes) {
				// System.out.println(chromosome.getSAMString());
				String chrName = chromosome.getSequenceName();
				int start = chromosome.getStart();
				int end = chromosome.getEnd();
				tasks.addAll(getChromosomeLevelTasks(records, chrName, start, end));
			}
		} catch(Exception e) {
			e.printStackTrace();
			System.exit(1);
		}
		
		// all tasks should be lib count.
		for(Task task : tasks) {
			task.type = Constants.TYPE_TARGET_MODE_LIBRARY_ESTIMATION_TASK;
		}
		
		return tasks;
	}

	@Override
	public int compareTo(Task o) {
		// TODO Auto-generated method stub
		if(this.type < o.type) {
			return 1;
		} else if(this.type > o.type) {
			return -1;
		} else if(this.readType < o.readType) {
			return 1;
		} else if(this.readType > o.readType) {
			return -1;
		}
		
		return 0;
	}
}
