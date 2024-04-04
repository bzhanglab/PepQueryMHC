package progistar.scan.function;

import java.io.File;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Hashtable;

import org.ahocorasick.trie.Emit;
import org.ahocorasick.trie.Trie;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import progistar.scan.data.BAMSRecord;
import progistar.scan.data.Constants;
import progistar.scan.run.Scan;
import progistar.scan.run.Task;

public class TargetModeRun {

	

	public static void runTargetMode (Task task) {
		if(task.type == Constants.TYPE_MAPPED_TASK) {
			countMappedReads(task);
		} else if(task.type == Constants.TYPE_UNMAPPED_TASK) {
			countUnmappedReads(task);
		}
	}


	/**
	 * If mode is all, than it uses all records.
	 * On the other hand, it uses only records with unmapped reads.
	 * 
	 * 
	 * @param records
	 * @param mode
	 * @return
	 */
	private static Trie getTrie (ArrayList<BAMSRecord> records) {
		ArrayList<String> sequences = new ArrayList<String>();
		Hashtable<String, String> rmDups = new Hashtable<String, String>();
		
		
		for(BAMSRecord record : records) {
			if(rmDups.get(record.sequence) == null) {
				sequences.add(record.sequence);
				rmDups.put(record.sequence, "");
			}
		}
		
		if(sequences.size() == 0) {
			return null;
		} else {
			return Trie.builder().addKeywords(sequences).build();
		}
	}
	

	private static void countUnmappedReads(Task task) {
		long startTime = System.currentTimeMillis();
		try (SamReader samReader = SamReaderFactory.makeDefault().open(new File(Scan.inputFilePath))) {
			// for unmapped reads
			Trie trie = getTrie(task.records);
			
			Hashtable<String, Integer> totalCounts = new Hashtable<String, Integer>();
			
			SAMRecordIterator iterator = samReader.queryUnmapped();
			while (iterator.hasNext()) {
                SAMRecord samRecord = iterator.next();
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
            		findings.put(emit.getKeyword(), "");
            	}
            	
            	// increase +1
            	findings.forEach((peptide, nil)->{
            		Integer cnt = totalCounts.get(peptide);
            		if(cnt == null) {
            			cnt = 0;
            		}
            		totalCounts.put(peptide, cnt+1)
;            	});
            	
            }
            iterator.close();
            
            task.records.forEach(record -> {
            	Integer cnt = totalCounts.get(record.sequence);
            	if(cnt == null) {
            		cnt = 0;
            	}
            	
            	record.readCnt += cnt;
            });
            
		} catch(Exception e) {
			e.printStackTrace();
			System.exit(1);
		}
		long endTime = System.currentTimeMillis();
		System.out.println((endTime-startTime)/1000+" sec");
	}
	
	
	private static void countMappedReads (Task task) {
		long startTime = System.currentTimeMillis();
		try (SamReader samReader = SamReaderFactory.makeDefault().open(new File(Scan.inputFilePath))) {
			double size = task.records.size();
			
			for(int i=0; i<size; i++) {

				BAMSRecord record = task.records.get(i);
				Trie trie = Trie.builder().addKeyword(record.sequence).build();
				
				SAMRecordIterator iterator = samReader.queryOverlapping(record.chr, record.start, record.end);
				int leastCount = find(iterator, record, trie);
	            
	            //
	            // in case of soft-clip, it can be zero because of unstable record range.
	            //////////////////////////////////// START SFOT-CLIP /////////////////////
	            if(leastCount == 0) {
	            	iterator = samReader.queryOverlapping(record.chr, record.start-100, record.end+100);
	            	find(iterator, record, trie);
	            }
	            ////////////////////////////////////// END SFOT-CLIP /////////////////////
			}
            
		} catch(Exception e) {
			e.printStackTrace();
			System.exit(1);
		}
		long endTime = System.currentTimeMillis();
		System.out.println((endTime-startTime)/1000+" sec");
	}
	
	private static int find (SAMRecordIterator iterator, BAMSRecord record, Trie trie) {
		int leastCount = 0;
		while (iterator.hasNext()) {
            SAMRecord samRecord = iterator.next();
            // Process each SAM record
            if(samRecord.getStart() <= record.start && samRecord.getEnd() >= record.end) {
                String sequence = samRecord.getReadString();
                boolean isFound = false;
                if(Scan.sequence.equalsIgnoreCase(Constants.SEQUENCE_NUCLEOTIDE)) {
                	if(trie.parseText(sequence).size() > 0) {
                		isFound = true;
                	}
                } else if(Scan.sequence.equalsIgnoreCase(Constants.SEQUENCE_PEPTIDE)) {
                	if(record.strand.equalsIgnoreCase("-")) {
                		sequence = Translator.getReverseComplement(sequence);
                	}
                	
                	for(int fr=0; fr<3; fr++) {
                		String peptide = Translator.translation(sequence, fr);
                		if(trie.parseText(peptide).size() > 0) {
	                		isFound = true;
	                		break;
	                	}
                	}
                }
                
                if(isFound) {
                	record.readCnt++;
                	leastCount++;
                }
           }
            
        }
        iterator.close();
        
        return leastCount;
	}
}
