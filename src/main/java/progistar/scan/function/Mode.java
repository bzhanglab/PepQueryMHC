package progistar.scan.function;

import java.util.ArrayList;
import java.util.Collection;

import org.ahocorasick.trie.Emit;
import org.ahocorasick.trie.Trie;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.fastq.FastqRecord;
import progistar.scan.data.BarcodeTable;
import progistar.scan.data.Constants;
import progistar.scan.data.LocTable;
import progistar.scan.data.LocationInformation;
import progistar.scan.data.Parameters;
import progistar.scan.run.Main;
import progistar.scan.run.Task;

public abstract class Mode {
	
	private static ArrayList<Character> getStrandedness (int flags) {
		ArrayList<Character> strands = new ArrayList<Character>();
		boolean isFirstSegment = (0x40 & flags) == 0x40 ? true : false;
		boolean isForward = (0x10 & flags) == 0x10 ? false : true;
		
		// non-stranded
		if(Parameters.strandedness.equalsIgnoreCase(Constants.NON_STRANDED)) {
			strands.add('+');
			strands.add('-');
		}  
		// Single-end
		else if(Parameters.strandedness.equalsIgnoreCase(Constants.F_STRANDED)) {
			if(isForward) {
				strands.add('+');
			} else {
				strands.add('-');
			}
		} else if(Parameters.strandedness.equalsIgnoreCase(Constants.R_STRANDED)) {
			if(isForward) {
				strands.add('-');
			} else {
				strands.add('+');
			}
		} 
		// Paired-end
		else {
			// R1
			if(isFirstSegment) {
				if(Parameters.strandedness.equalsIgnoreCase(Constants.FR_STRANDED)) {
					if(isForward) {
						strands.add('+');
					} else {
						strands.add('-');
					}
				} else if(Parameters.strandedness.equalsIgnoreCase(Constants.RF_STRANDED)) {
					if(isForward) {
						strands.add('-');
					} else {
						strands.add('+');
					}
				}
			} 
			// R2
			else {
				if(Parameters.strandedness.equalsIgnoreCase(Constants.FR_STRANDED)) {
					if(isForward) {
						strands.add('-');
					} else {
						strands.add('+');
					}
				} else if(Parameters.strandedness.equalsIgnoreCase(Constants.RF_STRANDED)) {
					if(isForward) {
						strands.add('+');
					} else {
						strands.add('-');
					}
				} 
			}
		}
		
		
		
		return strands;
	}
	
	/**
	 * For SCAN/TARGET modes.
	 * Task is an independent job, one thread takes only one task.
	 * Plus, task is not shared between threads.
	 * 
	 * @param iterator
	 * @param trie
	 * @param task
	 */
	public static void find (SAMRecordIterator iterator, Trie trie, Task task) {
		int count = 0;
		while (iterator.hasNext()) {
            SAMRecord samRecord = iterator.next();
            count ++;
            
            boolean isPass = false;
            // if the task is for mapped reads
            // only reads with below that genomic start are retrieved
            if(task.readType == Constants.MAPPED_READS) {
            	if(Parameters.count.equalsIgnoreCase(Constants.COUNT_PRIMARY) && samRecord.isSecondaryAlignment()) {
            		isPass = true;
            	}
            	
            	// In case of ScanMode, it should check the task range
            	// In case of TargetMode, iterator is already checked by a previous call.
            	if(task.type == Constants.TYPE_SCAN_MODE_TASK) {
            		if( !(samRecord.getAlignmentStart() >= task.start && 
            				samRecord.getAlignmentStart() <= task.end) ) {
            			isPass = true;
            		}
            	}
            }
            // if the task is for unmapped reads
            // only process reads given range
            else if(task.readType == Constants.UNMAPPED_READS) {
            	if(count < task.start || count > task.end) {
        			isPass = true;
        		}
            }

            
            if(isPass) {
            	continue;
            }
            
            // if barcode id is null or others, pass the read
            // no worry about bulk RNA-seq because it should be "undefined" in the bulk RNA-seq.
            String barcodeId = BarcodeTable.getBarcodeFromBam(samRecord);
            if(barcodeId.equalsIgnoreCase(Constants.NULL_BARCODE_ID) || barcodeId.equalsIgnoreCase(Constants.OTHER_BARCODE_ID) ) {
            	continue;
            }
            
            
            // increase processed reads if and only if
            // a read is primary.
            // In case of target mode, we do not count the reads in this routine.
            if(!samRecord.isSecondaryAlignment() && task.type == Constants.TYPE_SCAN_MODE_TASK) {
            	Double pReads = task.processedReads.get(barcodeId);
            	if(pReads == null) {
            		pReads = .0;
            	}
            	pReads++;
            	task.processedReads.put(barcodeId, pReads);
            }
            
            // Process each SAM record
            // determine strand
            int flags = samRecord.getFlags();
            ArrayList<Character> strands = getStrandedness(flags);
            
            for(Character strand : strands) {
            	String sequence = null;
            	if(strand == '+') {
            		sequence = samRecord.getReadString();
            	} else {
            		sequence = Translator.getReverseComplement(samRecord.getReadString());
            	}
            	
            	if(Parameters.sequence.equalsIgnoreCase(Constants.SEQUENCE_NUCLEOTIDE)) {
            		Collection<Emit> emits = trie.parseText(sequence);
            		
            		for(Emit emit : emits) {
        				LocationInformation matchedLocation = LocationInformation.getMatchedLocation(samRecord, emit, 0, strand);
        				if(matchedLocation != null) {
        					matchedLocation.inputSequence = emit.getKeyword();
        					if(task.type == Constants.TYPE_TARGET_MODE_TASK) {
        						
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
            	} else if(Parameters.sequence.equalsIgnoreCase(Constants.SEQUENCE_PEPTIDE)) {
            		for(int fr=0; fr<3; fr++) {
            			String peptide = Translator.translation(sequence, fr);
            			// if it is il equal mode?
            			if(Parameters.isILEqual) {
            				peptide = peptide.replace("I", "L");
            			}
            			
            			Collection<Emit> emits = trie.parseText(peptide);
            			
            			for(Emit emit : emits) {
            				LocationInformation matchedLocation = LocationInformation.getMatchedLocation(samRecord, emit, fr, strand);
            				if(matchedLocation != null) {
            					matchedLocation.inputSequence = emit.getKeyword();
            					
            					// we are only interested in the given region in case of target mode.
            					if(task.type == Constants.TYPE_TARGET_MODE_TASK) {
            						
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
	
	/**
	 * For FASTQ mode.
	 * Note that, task is a shared object between threads.
	 * Do not implement an operation using a task object that may cause concurrence between threads.
	 * 
	 * @param records
	 * @param trie
	 * @param task
	 * @return
	 */
	public static LocTable find (ArrayList<FastqRecord> records, Trie trie, Task task) {
		// Reads from a FASTQ file must be forwarded.
		int flags = 0x00;
		if(task.start == 0 || task.start == 1) {
			flags |= 0x40; // first segment
		} else if(task.start == 2) {
			flags |= 0x80; // last segment
		}
        ArrayList<Character> strands = getStrandedness(flags);
        LocTable locTable = new LocTable();
        
        for(FastqRecord fastqRecord : records) {
        	for(Character strand : strands) {
            	String sequence = null;
            	if(strand == '+') {
            		sequence = fastqRecord.getReadString();
            	} else {
            		sequence = Translator.getReverseComplement(fastqRecord.getReadString());
            	}
            	
            	if(Parameters.sequence.equalsIgnoreCase(Constants.SEQUENCE_NUCLEOTIDE)) {
            		Collection<Emit> emits = trie.parseText(sequence);
            		
            		for(Emit emit : emits) {
        				LocationInformation matchedLocation = LocationInformation.getMatchedLocation(fastqRecord, emit, 0, strand);
        				if(matchedLocation != null) {
        					matchedLocation.inputSequence = emit.getKeyword();
        					if(locTable.putLocation(matchedLocation)) {
        						matchedLocation.calMetaInfo();
        					}
        				}
        			}
            	} else if(Parameters.sequence.equalsIgnoreCase(Constants.SEQUENCE_PEPTIDE)) {
            		for(int fr=0; fr<3; fr++) {
            			String peptide = Translator.translation(sequence, fr);
            			// if it is il equal mode?
            			if(Parameters.isILEqual) {
            				peptide = peptide.replace("I", "L");
            			}
            			
            			Collection<Emit> emits = trie.parseText(peptide);
            			
            			for(Emit emit : emits) {
            				LocationInformation matchedLocation = LocationInformation.getMatchedLocation(fastqRecord, emit, fr, strand);
            				if(matchedLocation != null) {
            					matchedLocation.inputSequence = emit.getKeyword();
            					if(locTable.putLocation(matchedLocation)) {
            						matchedLocation.calMetaInfo();
            					}
            				}
            			}
            		}
            	}
            }
        }
        
        return locTable;
	}
}
