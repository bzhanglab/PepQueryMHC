package progistar.scan.function;

import java.io.File;
import java.util.Collection;

import org.ahocorasick.trie.Emit;
import org.ahocorasick.trie.Trie;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import progistar.scan.data.Constants;
import progistar.scan.data.LocationInformation;
import progistar.scan.run.Scan;
import progistar.scan.run.Task;

public class ScanModeRun {
	
	public static void runScanMode (Task task) {
		if(task.type == Constants.TYPE_SCAN_MODE_TASK) {
			if(Scan.verbose) {
				System.out.println(task.chrName+":"+task.start+"-"+task.end);
			}
			scanReads(task);
		}
	}
	
	private static void scanReads (Task task) {
		long startTime = System.currentTimeMillis();
		// to prevent racing
		File file = new File(Scan.bamFile.getAbsolutePath());
		try (SamReader samReader = SamReaderFactory.makeDefault().open(file)) {
			SAMRecordIterator iterator = null;
			if(task.readType == Constants.MAPPED_READS) {
				iterator = samReader.query(task.chrName, task.start, task.end, false);
			} else {
				iterator = samReader.queryUnmapped();
			}
			find(iterator, Task.allTrie, task);
		} catch(Exception e) {
			e.printStackTrace();
			System.exit(1);
		}
		
		long endTime = System.currentTimeMillis();
		if(Scan.verbose) {
			System.out.println("Task"+task.taskIdx+" "+(endTime-startTime)/1000+" sec");
		}
	}
	
	
	protected static void find (SAMRecordIterator iterator, Trie trie, Task task) {
		int count = 0;
		while (iterator.hasNext()) {
            SAMRecord samRecord = iterator.next();
            count ++;
            
            boolean isPass = false;
            // if the task is for mapped reads
            // only reads with below that genomic start are retrieved
            if(task.readType == Constants.MAPPED_READS) {
            	if(Scan.count.equalsIgnoreCase(Constants.COUNT_PRIMARY) && samRecord.isSecondaryAlignment()) {
            		isPass = true;
            	}
            	
            	// In case of ScanMode, it should check the task range
            	// In case of TargetMode, iterator is already checked by a previous call.
            	if(task.type == Constants.TYPE_SCAN_MODE_TASK) {
            		if( !(samRecord.getAlignmentStart() >= task.start && 
            				samRecord.getAlignmentStart() < task.end) ) {
            			isPass = true;
            		}
            	}
            	
            	
            }
            // if the task is for unmapped reads
            // only process reads given range
            else {
            	// TODO: validate this logic in target mode
            	if(task.type == Constants.TYPE_SCAN_MODE_TASK) {
            		if(count < task.start || count > task.end) {
            			isPass = true;
            		}
            	}
            }

            
            if(isPass) {
            	continue;
            }
            
            // increase processed reads if and only if
            // a read is primary and mapped to a genomic region.
            if(!samRecord.isSecondaryAlignment()) {
            	task.processedReads++;
            }
            
            // Process each SAM record
            String[] frrvSequences = {samRecord.getReadString(), Translator.getReverseComplement(samRecord.getReadString())};
            
            int strand = 0;
            for(String sequence : frrvSequences) {
            	char strandChar = '+';
            	if(strand != 0) {
            		strandChar = '-';
            	}
            	strand++;
            	
            	if(Scan.sequence.equalsIgnoreCase(Constants.SEQUENCE_NUCLEOTIDE)) {
            		Collection<Emit> emits = trie.parseText(sequence);
            		
            		for(Emit emit : emits) {
        				LocationInformation matchedLocation = LocationInformation.getMatchedLocation(samRecord, emit, 0, strandChar);
        				if(matchedLocation != null) {
        					matchedLocation.inputSequence = emit.getKeyword();
        					if(task.type == Constants.TYPE_TARGET_MODE_MAPPED_TASK) {
        						
        						// discard if the location is not matched
        						int targetMappedTaskCurrentIdx = task.currentRecordIdx;
        						if(!task.records.get(targetMappedTaskCurrentIdx).location
        								.equalsIgnoreCase(matchedLocation.location)) {
        							continue;
        						}
        						
        					}
        					
        					if(task.locTable.putLocation(matchedLocation)) {
        						matchedLocation.calMetaInfo();
        					}
        				}
        			}
            	} else if(Scan.sequence.equalsIgnoreCase(Constants.SEQUENCE_PEPTIDE)) {
            		for(int fr=0; fr<3; fr++) {
            			String peptide = Translator.translation(sequence, fr);
            			// if it is il equal mode?
            			if(Scan.isILEqual) {
            				peptide = peptide.replace("I", "L");
            			}
            			
            			Collection<Emit> emits = trie.parseText(peptide);
            			
            			for(Emit emit : emits) {
            				LocationInformation matchedLocation = LocationInformation.getMatchedLocation(samRecord, emit, fr, strandChar);
            				if(matchedLocation != null) {
            					matchedLocation.inputSequence = emit.getKeyword();
            					if(task.type == Constants.TYPE_TARGET_MODE_MAPPED_TASK) {
            						
            						// discard if the location is not matched
            						int targetMappedTaskCurrentIdx = task.currentRecordIdx;
            						if(!task.records.get(targetMappedTaskCurrentIdx).location
            								.equalsIgnoreCase(matchedLocation.location)) {
            							continue;
            						}
            						
            					}
            					
            					if(task.locTable.putLocation(matchedLocation)) {
            						matchedLocation.calMetaInfo();
            					}
            				}
            			}
            		}
            	}
            }
            
        }
        iterator.close();
	}
	
}











