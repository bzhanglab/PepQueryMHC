package progistar.scan.run;

import java.io.File;
import java.util.ArrayList;
import java.util.Hashtable;
import java.util.List;

import org.ahocorasick.trie.Trie;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import progistar.scan.data.Constants;
import progistar.scan.data.LocTable;
import progistar.scan.data.Parameters;
import progistar.scan.data.SequenceRecord;

public class Task implements Comparable<Task> {

	public int taskIdx = -1;
	public int type = Constants.TYPE_TARGET_MODE_TASK;
	public int readType = Constants.MAPPED_READS;
	public Hashtable<String, Double> processedReads = new Hashtable<String, Double>();
	// only available for TargetMode
	public ArrayList<SequenceRecord> records = new ArrayList<SequenceRecord>();
	public int currentRecordIdx = 0;
	
	// only available for ScanMode.
	public static Trie allTrie = null;
	public LocTable locTable = new LocTable();
	public String chrName;
	public int start;
	public int end;
	
	
	// estimated memory usage
	public long peakMemory = 0;
	
	// for strand type
	public int R1F = 0;
	public int R1R = 0;
	public int R2F = 0;
	public int R2R = 0;
	
	public String getTaskInfo () {
		String typeStr = null;
		if(type == Constants.TYPE_TARGET_MODE_TASK) {
			if(readType == Constants.MAPPED_READS) {
				typeStr = "target mode (mapped)";
			} else {
				typeStr = "target mode (unmapped)";
			}
		} else if(type == Constants.TYPE_SCAN_MODE_TASK) {
			typeStr = "scan mode";
		} else if(type == Constants.TYPE_TARGET_MODE_LIBRARY_ESTIMATION_TASK) {
			typeStr = "target mode (libary estimation)";
		} else if(type == Constants.TYPE_STRAND_DETECTION_TASK) {
			typeStr = "strand detection";
		} else if(type == Constants.TYPE_EXTRACT_MODE_TASK) {
			typeStr = "extract mode";
		}
		
		return typeStr+": Task"+taskIdx;
	}
	
	public Task(int type) {
		this.type = type;
	}
	
	private static ArrayList<Task> getStrandDectectionTasks (String chrName, int size, int mode) {
		assert mode == Constants.TYPE_STRAND_DETECTION_TASK;
		
		ArrayList<Task> tasks = new ArrayList<Task>();
		Task task = new Task(mode);
		task.readType = Constants.MAPPED_READS;
		task.records = null;
		task.chrName = chrName;
		task.start = 1;
		task.end = size;
		tasks.add(task);
		
		return tasks;
	}
	
	private static ArrayList<Task> getChromosomeLevelTasks (ArrayList<SequenceRecord> records, 
			String chrName, int start, int end, int mode) {
		assert mode == Constants.TYPE_SCAN_MODE_TASK || mode == Constants.TYPE_TARGET_MODE_TASK;
		
		ArrayList<Task> tasks = new ArrayList<Task>();
		
		int divider = Parameters.threadNum * Parameters.chunkSize;
		// if unmapped reads?
		// only number of threads will be partitioned
		if(chrName.equalsIgnoreCase(Constants.NULL)) {
			divider = Parameters.threadNum;
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
			
			Task task = new Task(mode);
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
		
		File file = new File(Parameters.bamFile.getAbsolutePath());
		// build global trie
		System.out.println("Build Trie");
		Task.allTrie = SequenceRecord.getTrie(records);
		System.out.println("Complete building Trie");
		try (SamReader samReader = SamReaderFactory.makeDefault().open(file)) {
			List<SAMSequenceRecord> chromosomes = samReader.getFileHeader().getSequenceDictionary().getSequences();
			for(SAMSequenceRecord chromosome : chromosomes) {
				// System.out.println(chromosome.getSAMString());
				String chrName = chromosome.getSequenceName();
				int start = chromosome.getStart();
				int end = chromosome.getEnd();
				
				tasks.addAll(getChromosomeLevelTasks(records, chrName, start, end, Constants.TYPE_SCAN_MODE_TASK));
			}
			
			// for unmapped reads
			SAMRecordIterator unmappedIter = samReader.queryUnmapped();
			int size = 0;
			while(unmappedIter.hasNext()) {
				SAMRecord samRecord = unmappedIter.next();
				if(Parameters.unmmapedMarker == null) {
					Parameters.unmmapedMarker = samRecord.getReferenceName();
				}
				size ++;
			}
			if(Parameters.unmmapedMarker != null) {
				// System.out.println("@SQ\t"+Scan.unmmapedMarker+"\tLN:"+size);
				tasks.addAll(getChromosomeLevelTasks(records, Constants.NULL, 1, size, Constants.TYPE_SCAN_MODE_TASK));
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
		System.out.println("Unmapped records: "+unmappedSize);
		
		// generate unmapped tasks
		File file = new File(Parameters.bamFile.getAbsolutePath());
		if(unmappedSize != 0) {
			try (SamReader samReader = SamReaderFactory.makeDefault().open(file)) {
				// for unmapped reads
				SAMRecordIterator unmappedIter = samReader.queryUnmapped();
				int size = 0;
				while(unmappedIter.hasNext()) {
					SAMRecord samRecord = unmappedIter.next();
					if(Parameters.unmmapedMarker == null) {
						Parameters.unmmapedMarker = samRecord.getReferenceName();
					}
					size ++;
				}
				if(Parameters.unmmapedMarker != null) {
					// System.out.println("@SQ\t"+Scan.unmmapedMarker+"\tLN:"+size);
					tasks.addAll(getChromosomeLevelTasks(unmappedRecords, Constants.NULL, 1, size, Constants.TYPE_TARGET_MODE_TASK));
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
		}
		
		
		// generate mapped tasks
		int sIdx = 0;
		while(sIdx < mappedSize) {
			int eIdx = sIdx + chunkSize > mappedSize ? mappedSize : sIdx + chunkSize;
			Task task = new Task(Constants.TYPE_TARGET_MODE_TASK);
			task.readType = Constants.MAPPED_READS;
			
			for(int i=sIdx; i<eIdx; i++) {
				task.records.add(mappedRecords.get(i));
			}
			tasks.add(task);
			task.taskIdx = tasks.size();
			sIdx = eIdx;
		}
		
		if(Parameters.verbose) {
			System.out.println("Task list");
			for(Task task : tasks) {
				System.out.println(task.getTaskInfo());
			}
		}
		
		return tasks;
	}
	
	public static ArrayList<Task> getFASTQModeTasks (ArrayList<SequenceRecord> records) {
		ArrayList<Task> tasks = new ArrayList<Task>();
		
		// build global trie
		System.out.println("Build Trie");
		Task.allTrie = SequenceRecord.getTrie(records);
		System.out.println("Complete building Trie");
			
		
		if(Parameters.sequencingFileType == Constants.SEQ_FASTQ_SINGLE) {
			Task task = new Task(Constants.TYPE_FASTQ_MODE_TASK);
			task.start = 0;
			tasks.add(task);
		} else if(Parameters.sequencingFileType == Constants.SEQ_FASTQ_PAIRED) {
			// first fastq
			Task task = new Task(Constants.TYPE_FASTQ_MODE_TASK);
			task.start = 1;
			tasks.add(task);
			// second fastq
			task = new Task(Constants.TYPE_FASTQ_MODE_TASK);
			task.start = 2;
			tasks.add(task);
		}
		
		// assign idx
		for(int i=0; i<tasks.size(); i++) {
			tasks.get(i).taskIdx = (i+1);
		}
		
		return tasks;
	}
	
	public static ArrayList<Task> getStrandDetectionTask () {
		ArrayList<Task> tasks = new ArrayList<Task>();
		
		File file = new File(Parameters.bamFile.getAbsolutePath());
		Task.allTrie = null;
		
		try (SamReader samReader = SamReaderFactory.makeDefault().open(file)) {
			// System.out.println(samReader.getFileHeader().getSequenceDictionary().getSequences().get(0).getSequenceLength());
			List<SAMSequenceRecord> chromosomes = samReader.getFileHeader().getSequenceDictionary().getSequences();
			for(SAMSequenceRecord chromosome : chromosomes) {
				// System.out.println(chromosome.getSAMString());
				String chrName = chromosome.getSequenceName();
				tasks.addAll(getStrandDectectionTasks(chrName, 100000, Constants.TYPE_STRAND_DETECTION_TASK));
			}
			
			// assign idx
			for(int i=0; i<tasks.size(); i++) {
				tasks.get(i).taskIdx = (i+1);
			}
		} catch(Exception e) {
			e.printStackTrace();
			System.exit(1);
		}
		
		return tasks;
	}
	
	public static ArrayList<Task> getLibSizeTask () {
		ArrayList<Task> tasks = new ArrayList<Task>();
		
		File file = new File(Parameters.bamFile.getAbsolutePath());
		Task.allTrie = null;
		
		try (SamReader samReader = SamReaderFactory.makeDefault().open(file)) {
			// System.out.println(samReader.getFileHeader().getSequenceDictionary().getSequences().get(0).getSequenceLength());
			List<SAMSequenceRecord> chromosomes = samReader.getFileHeader().getSequenceDictionary().getSequences();
			for(SAMSequenceRecord chromosome : chromosomes) {
				// System.out.println(chromosome.getSAMString());
				String chrName = chromosome.getSequenceName();
				int start = chromosome.getStart();
				int end = chromosome.getEnd();
				tasks.addAll(getChromosomeLevelTasks(null, chrName, start, end, Constants.TYPE_TARGET_MODE_TASK));
			}
			
			// for unmapped reads
			SAMRecordIterator unmappedIter = samReader.queryUnmapped();
			int size = 0;
			while(unmappedIter.hasNext()) {
				SAMRecord samRecord = unmappedIter.next();
				if(Parameters.unmmapedMarker == null) {
					Parameters.unmmapedMarker = samRecord.getReferenceName();
				}
				size ++;
			}
			if(Parameters.unmmapedMarker != null) {
				// System.out.println("@SQ\t"+Scan.unmmapedMarker+"\tLN:"+size);
				tasks.addAll(getChromosomeLevelTasks(null, Constants.NULL, 1, size, Constants.TYPE_TARGET_MODE_TASK));
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
	
	public static ArrayList<Task> getExtractModeTasks (ArrayList<SequenceRecord> records, int chunkSize) {
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
		System.out.println("Unmapped records: "+unmappedSize);
		
		// generate unmapped tasks
		File file = new File(Parameters.bamFile.getAbsolutePath());
		if(unmappedSize != 0) {
			try (SamReader samReader = SamReaderFactory.makeDefault().open(file)) {
				// for unmapped reads
				SAMRecordIterator unmappedIter = samReader.queryUnmapped();
				int size = 0;
				while(unmappedIter.hasNext()) {
					SAMRecord samRecord = unmappedIter.next();
					if(Parameters.unmmapedMarker == null) {
						Parameters.unmmapedMarker = samRecord.getReferenceName();
					}
					size ++;
				}
				if(Parameters.unmmapedMarker != null) {
					// System.out.println("@SQ\t"+Scan.unmmapedMarker+"\tLN:"+size);
					tasks.addAll(getChromosomeLevelTasks(unmappedRecords, Constants.NULL, 1, size, Constants.TYPE_EXTRACT_MODE_TASK));
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
		}
		
		
		// generate mapped tasks
		int sIdx = 0;
		while(sIdx < mappedSize) {
			int eIdx = sIdx + chunkSize > mappedSize ? mappedSize : sIdx + chunkSize;
			Task task = new Task(Constants.TYPE_EXTRACT_MODE_TASK);
			task.readType = Constants.MAPPED_READS;
			
			for(int i=sIdx; i<eIdx; i++) {
				task.records.add(mappedRecords.get(i));
			}
			tasks.add(task);
			task.taskIdx = tasks.size();
			sIdx = eIdx;
		}
		
		if(Parameters.verbose) {
			System.out.println("Task list");
			for(Task task : tasks) {
				System.out.println(task.getTaskInfo());
			}
		}
		
		return tasks;
	}
	
	public String getExtractedBamFileNameForUnmapped () {
		return "u.tmp."+this.taskIdx+".bam";
	}
	
	public String getExtractedBamFileNameForMapped () {
		return "m.tmp."+this.taskIdx+".bam";
	}
}
