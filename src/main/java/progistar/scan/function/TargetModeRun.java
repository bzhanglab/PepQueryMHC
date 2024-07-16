package progistar.scan.function;

import java.io.File;
import java.util.Collection;
import java.util.Hashtable;

import org.ahocorasick.trie.Emit;
import org.ahocorasick.trie.Trie;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import progistar.scan.data.BarcodeTable;
import progistar.scan.data.Constants;
import progistar.scan.data.LocationInformation;
import progistar.scan.data.SequenceRecord;
import progistar.scan.run.Scan;
import progistar.scan.run.Task;

public class TargetModeRun {
	
	public static void runTargetMode (Task task) {
		if(task.type == Constants.TYPE_TARGET_MODE_MAPPED_TASK) {
			countMappedReads(task);
		} else if(task.type == Constants.TYPE_TARGET_MODE_UNMAPPED_TASK) {
			countUnmappedReads(task);
		} else if(task.type == Constants.TYPE_TARGET_MODE_LIBRARY_ESTIMATION_TASK) {
			estimateLibSize(task);
		}
	}

	private static void countUnmappedReads(Task task) {
		long startTime = System.currentTimeMillis();
		// to prevent racing
		File file = new File(Scan.bamFile.getAbsolutePath());
		try (SamReader samReader = SamReaderFactory.makeDefault().open(file)) {
			// for unmapped reads
			Trie trie = SequenceRecord.getTrie(task.records);
			
			Hashtable<String, Hashtable<String, Long>> totalCounts = new Hashtable<String, Hashtable<String, Long>>();
			SAMRecordIterator iterator = samReader.queryUnmapped();
			while (iterator.hasNext()) {
                SAMRecord samRecord = iterator.next();
                
                if(Scan.count.equalsIgnoreCase(Constants.COUNT_PRIMARY) && samRecord.isSecondaryAlignment()) {
                	continue;
                }
                
                String barcodeId = BarcodeTable.getBarcodeFromBam(samRecord);
                
                // Process each SAM record
                String frSequence = samRecord.getReadString();
                String rcSequence = Translator.getReverseComplement(frSequence);
                Hashtable<String, String> findings = new Hashtable<String, String>();
                Collection<Emit> emits = null;
                if(Scan.sequence.equalsIgnoreCase(Constants.SEQUENCE_NUCLEOTIDE)) {
                	emits = trie.parseText(frSequence);
                	emits.addAll(trie.parseText(rcSequence));
                } else if(Scan.sequence.equalsIgnoreCase(Constants.SEQUENCE_PEPTIDE)) {
                	for(int fr=0; fr<3; fr++) {
                		String peptide = Translator.translation(frSequence, fr);
                		if(emits == null) {
                			emits = trie.parseText(peptide);
                		} else {
                			emits.addAll(trie.parseText(peptide));
                		}
                		
                		peptide = Translator.translation(rcSequence, fr);
                		emits.addAll(trie.parseText(peptide));
                	}
                }
                
                // save findings
            	for(Emit emit : emits) {
            		findings.put(emit.getKeyword(), barcodeId);
            	}
            	
            	// increase +1
            	findings.forEach((peptide, id)->{
            		Hashtable<String, Long> counts = totalCounts.get(peptide);
            		if(counts == null) {
            			counts = new Hashtable<String, Long>();
            			totalCounts.put(peptide, counts);
            		}
            		
            		Long val = counts.get(id);
            		if(val == null) {
            			val = 0L;
            		}
            		counts.put(id, val + 1);
            	});
            	
            }
            iterator.close();
            
            task.records.forEach(record -> {
            	Hashtable<String, Long> counts = totalCounts.get(record.sequence);
            	if(counts != null) {
            		record.readCounts = counts;
            	}
            });
            
		} catch(Exception e) {
			e.printStackTrace();
			System.exit(1);
		}
		long endTime = System.currentTimeMillis();
		
		if(Scan.verbose) {
			System.out.println("Task"+task.taskIdx+" "+(endTime-startTime)/1000+" sec");
		}
	}
	
	
	private static void countMappedReads (Task task) {
		long startTime = System.currentTimeMillis();
		// to prevent racing
		File file = new File(Scan.bamFile.getAbsolutePath());
		try (SamReader samReader = SamReaderFactory.makeDefault().open(file)) {
			double size = task.records.size();
			for(int i=0; i<size; i++) {

				SequenceRecord record = task.records.get(i);
				Trie trie = Trie.builder().addKeyword(record.sequence).build();
				
				SAMRecordIterator iterator = samReader.queryOverlapping(record.chr, record.start, record.end);
				int leastCount = find(iterator, record, trie, true);
	            
	            //
	            // in case of soft-clip, it can be zero because of unstable record range.
	            //////////////////////////////////// START SFOT-CLIP /////////////////////
	            if(leastCount == 0) {
	            	iterator = samReader.queryOverlapping(record.chr, record.start-100, record.end+100);
	            	find(iterator, record, trie, false);
	            }
	            ////////////////////////////////////// END SFOT-CLIP /////////////////////
			}
            
		} catch(Exception e) {
			e.printStackTrace();
			System.exit(1);
		}
		long endTime = System.currentTimeMillis();
		
		if(Scan.verbose) {
			System.out.println(task.taskIdx+" "+(endTime-startTime)/1000+" sec");
		}
	}
	
	/**
	 * trie has a single sequence.
	 * 
	 * @param iterator
	 * @param record
	 * @param trie
	 * @param included
	 * @return
	 */
	private static int find (SAMRecordIterator iterator, SequenceRecord record, Trie trie, boolean included) {
		int leastCount = 0;
		char strand = record.strand.charAt(0);
		while (iterator.hasNext()) {
            SAMRecord samRecord = iterator.next();
            LocationInformation matchedLocation = null;
            
            
            if(Scan.count.equalsIgnoreCase(Constants.COUNT_PRIMARY) && samRecord.isSecondaryAlignment()) {
            	continue;
            }
            
        	String sequence = samRecord.getReadString();
            boolean isFound = false;
            if(Scan.sequence.equalsIgnoreCase(Constants.SEQUENCE_NUCLEOTIDE)) {
            	Collection<Emit> emits = trie.parseText(sequence);
        		
        		for(Emit emit : emits) {
    				matchedLocation = LocationInformation.getMatchedLocation(samRecord, emit, 0, strand);
    				matchedLocation.inputSequence = emit.getKeyword();
    				if(record.location.equalsIgnoreCase(matchedLocation.location)) {
    					isFound = true;
    					break;
    				}
    			}
            } else if(Scan.sequence.equalsIgnoreCase(Constants.SEQUENCE_PEPTIDE)) {
            	if(record.strand.equalsIgnoreCase("-")) {
            		sequence = Translator.getReverseComplement(sequence);
            	}
            	
            	for(int fr=0; fr<3; fr++) {
            		String peptide = Translator.translation(sequence, fr);
            		Collection<Emit> emits = trie.parseText(peptide);
            		
            		for(Emit emit : emits) {
        				matchedLocation = LocationInformation.getMatchedLocation(samRecord, emit, fr, strand);
        				if(record.location.equalsIgnoreCase(matchedLocation.location)) {
        					isFound = true;
        					fr = 3;
        					break;
        				}
        			}
            	}
            }
            
            if(isFound) {
            	matchedLocation.readCounts.forEach((barcodeId, value)->{
            		Long val = record.readCounts.get(barcodeId);
            		if(val == null) {
            			val = 0L;
            		}
            		record.readCounts.put(barcodeId, value + val);
            	});
            	leastCount++;
            }
            
            
            
        }
        iterator.close();
        
        return leastCount;
	}
	
	private static void estimateLibSize (Task task) {
		long startTime = System.currentTimeMillis();
		// to prevent racing
		File file = new File(Scan.bamFile.getAbsolutePath());
		try (SamReader samReader = SamReaderFactory.makeDefault().open(file)) {
			SAMRecordIterator iterator = null;
			if(task.readType == Constants.MAPPED_READS) {
				iterator = samReader.query(task.chrName, task.start, task.end, false);
				estimate(iterator, task);
			}
		} catch(Exception e) {
			e.printStackTrace();
			System.exit(1);
		}
		
		long endTime = System.currentTimeMillis();
		
		if(Scan.verbose) {
			System.out.println("Task"+task.taskIdx+" "+(endTime-startTime)/1000+" sec");
		}
	}
	

	private static void estimate (SAMRecordIterator iterator, Task task) {
		while (iterator.hasNext()) {
            SAMRecord samRecord = iterator.next();
            boolean isPass = false;

            // if the task is for mapped reads
            // only reads with below that genomic start are retrieved
            if(task.readType == Constants.MAPPED_READS) {
            	if( !(samRecord.getAlignmentStart() >= task.start && 
            			samRecord.getAlignmentStart() < task.end) ) {
            		isPass = true;
            	}
            	
            	if(samRecord.isSecondaryAlignment()) {
            		isPass = true;
            	}
            	
            } else {
            	isPass = true;
            }
            
            if(isPass) {
            	continue;
            }
            
            task.processedReads++;
            
        }
        iterator.close();
	}
}
