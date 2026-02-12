package progistar.scan.function;

import java.io.File;
import java.util.ArrayList;

import org.ahocorasick.trie.Trie;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import progistar.scan.data.Constants;
import progistar.scan.data.Parameters;
import progistar.scan.data.SequenceRecord;
import progistar.scan.run.Task;

public class ExtractModeRun extends Mode {
	
	public static void runExtractMode (Task task) {
		if(task.type == Constants.TYPE_EXTRACT_MODE_TASK) {
			if(task.readType == Constants.MAPPED_READS) {
				countMappedReads(task);
			} else if(task.readType == Constants.UNMAPPED_READS) {
				countUnmappedReads(task);
			}
		}
	}
	
	private static void countUnmappedReads(Task task) {
		ArrayList<SAMRecord> unmappedSAMRecords = new ArrayList<SAMRecord>();
		
		long startTime = System.currentTimeMillis();
		// to prevent racing
		File file = new File(Parameters.bamFile.getAbsolutePath());
		try (SamReader samReader = SamReaderFactory.makeDefault().open(file)) {
			// for unmapped reads
			Trie trie = SequenceRecord.getTrie(task.records);
			SAMRecordIterator iterator = samReader.queryUnmapped();
			unmappedSAMRecords.addAll(find(iterator, trie, task, true));
			
			
			// write an unmapped file
			SAMFileHeader header = samReader.getFileHeader();
			if(unmappedSAMRecords.size() > 0) {
				try {
					File extractBam = new File(task.getExtractedBamFileNameForUnmapped());
					SAMFileWriter writer = new SAMFileWriterFactory().makeBAMWriter(header, false, extractBam);
					
					for(SAMRecord record : unmappedSAMRecords) {
						writer.addAlignment(record);
					}
					
					writer.close();
				}catch(Exception e) {
					
				}
			}
		} catch(Exception e) {
			e.printStackTrace();
			System.exit(1);
		}
		long endTime = System.currentTimeMillis();
		
		if(Parameters.verbose) {
			System.out.println("Task"+task.taskIdx+" "+(endTime-startTime)/1000+" sec");
		}
		
		
	}
	
	private static void countMappedReads (Task task) {
		ArrayList<SAMRecord> mappedSAMRecords = new ArrayList<SAMRecord>();
		
		long startTime = System.currentTimeMillis();
		// to prevent racing
		File file = new File(Parameters.bamFile.getAbsolutePath());
		try (SamReader samReader = SamReaderFactory.makeDefault().open(file)) {
			double size = task.records.size();
			for(int i=0; i<size; i++) {
				task.currentRecordIdx = i;
				SequenceRecord record = task.records.get(i);
				
				Trie trie = Trie.builder().addKeyword(record.sequence).build();
				
	            // in case of soft-clip, it can be zero because of unstable record range.
				SAMRecordIterator iterator = samReader.queryOverlapping(record.chr, record.start-100, record.end+100);
				mappedSAMRecords.addAll(find(iterator, trie, task, true));
			}
            
			
			// write a mapped file
			SAMFileHeader header = samReader.getFileHeader();
			if(mappedSAMRecords.size() > 0) {
				try {
					File extractBam = new File(task.getExtractedBamFileNameForMapped());
					SAMFileWriter writer = new SAMFileWriterFactory().makeBAMWriter(header, false, extractBam);
					
					for(SAMRecord record : mappedSAMRecords) {
						writer.addAlignment(record);
					}
					
					writer.close();
				}catch(Exception e) {
					
				}
			}
		} catch(Exception e) {
			e.printStackTrace();
			System.exit(1);
		}
		long endTime = System.currentTimeMillis();
		
		if(Parameters.verbose) {
			System.out.println(task.taskIdx+" "+(endTime-startTime)/1000+" sec");
		}
		
	}
}
