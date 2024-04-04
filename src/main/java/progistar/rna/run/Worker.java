package progistar.rna.run;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Hashtable;

import org.ahocorasick.trie.Emit;
import org.ahocorasick.trie.Trie;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.util.Locatable;
import progistar.rna.data.Constants;
import progistar.rna.data.Record;

public class Worker extends Thread {

	private Task task;
	private int workerID = -1;
	
	public Worker (int workerID, Task task) {
		super();
		this.task = task;
		this.workerID = workerID;
	}
	
	public void run() {
		System.out.println(task.fileInfo.file.getName()+" is running by "+this.workerID);
		if(Scan.mode.equalsIgnoreCase(Constants.MODE_SCAN)) {
			doMatchSam(task);
		} else if(Scan.mode.equalsIgnoreCase(Constants.MODE_TARGET)) {
			doMatchSamWithGenomicLoci(task);
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
	public Trie getTrie (ArrayList<Record> records, String mode) {
		ArrayList<String> nucleotides = new ArrayList<String>();
		Hashtable<String, String> rmDups = new Hashtable<String, String>();
		
		
		for(Record record : records) {
			
			if(mode.equalsIgnoreCase(Constants.MODE_SCAN)) {
				if(rmDups.get(record.frSequence) == null) {
					nucleotides.add(record.frSequence);
					rmDups.put(record.frSequence, "");
				}
				if(rmDups.get(record.rcSequence) == null) {
					nucleotides.add(record.rcSequence);
					rmDups.put(record.rcSequence, "");
				}
			} else if(mode.equalsIgnoreCase(Constants.MODE_TARGET)) {
				if(record.chr.equalsIgnoreCase("*")) {
					if(rmDups.get(record.frSequence) == null) {
						nucleotides.add(record.frSequence);
						rmDups.put(record.frSequence, "");
					}
					if(rmDups.get(record.rcSequence) == null) {
						nucleotides.add(record.rcSequence);
						rmDups.put(record.rcSequence, "");
					}
				}
			}
			
			
		}
		if(mode.equalsIgnoreCase(Constants.MODE_SCAN)) {
			System.out.println("A total of "+rmDups.size() +" in the mapped/unmapped reads");
		} else {
			System.out.println("A total of "+rmDups.size() +" in the unampped reads");
		}
		
		if(nucleotides.size() == 0) {
			return null;
		} else {
			return Trie.builder().addKeywords(nucleotides).build();
		}
	}
	
	public Hashtable<String, ArrayList<Record>> getMapBySequence (ArrayList<Record> records) {
		Hashtable<String, ArrayList<Record>> map = new Hashtable<String, ArrayList<Record>>();
		
		for(Record record : records) {
			String[] seqs = {record.frSequence, record.rcSequence};
			
			for(String seq : seqs) {
				ArrayList<Record> matchedRecords = map.get(seq);
				if(matchedRecords == null) {
					matchedRecords = new ArrayList<Record>();
					map.put(seq, matchedRecords);
				}
				matchedRecords.add(record);
			}
		}
		
		return map;
	}
	
	public void doMatchFastq(Task task) {
		System.out.println("Currently, this function is not implemented.");
		System.exit(1);
	}
	
	public void doMatchSam(Task task) {
		long startTime = System.currentTimeMillis();
		Trie trie = getTrie(task.records, Scan.mode);
		Hashtable<String, ArrayList<Record>> map = getMapBySequence(task.records);
		try (SamReader samReader = SamReaderFactory.makeDefault().open(task.fileInfo.file)) {
			
			for(SAMRecord samRecord : samReader) {
				// count the number of reads
				String sequence = samRecord.getReadString();
				Collection<Emit> emits = trie.parseText(sequence);
				for(Emit emit : emits) {
					String matchedSequence = emit.getKeyword();
					ArrayList<Record> mRecords = map.get(matchedSequence);
					
					for(Record record : mRecords) {
						if(record.frSequence.equalsIgnoreCase(matchedSequence)) {
							record.frCnt++;
						} else if(record.rcSequence.equalsIgnoreCase(matchedSequence)) {
							record.rcCnt++;
						} else {
							System.out.println("Out of scope!");
						}
					}
				}
			}
			
		} catch(Exception e) {
			
		}
		long endTime = System.currentTimeMillis();
		System.out.println((endTime-startTime)/1000+" sec");
	}
	
	
	public void doMatchSamWithGenomicLoci(Task task) {
		long startTime = System.currentTimeMillis();
		String[] chr = new String[1];
		try (SamReader samReader = SamReaderFactory.makeDefault().open(task.fileInfo.file)) {
			// for mapped reads
			Hashtable<String, ArrayList<Record>> uniqueRecords = new Hashtable<String, ArrayList<Record>>();
			for(Record record : task.records) {
				// pass unmapped reads
				if(record.chr.equalsIgnoreCase("*")) continue;
				
				String key = record.chr+":"+record.start+"-"+record.end+"_"+record.frSequence;
				ArrayList<Record> dupRecords = uniqueRecords.get(key);
				if(dupRecords == null) {
					dupRecords = new ArrayList<Record>();
					uniqueRecords.put(key, dupRecords);
				}
				dupRecords.add(record);
			}
			
			double size = uniqueRecords.size();
			double[] progress = new double[2];
			progress[1] = 0.1;
			uniqueRecords.forEach((key, dupRecords)->{
				progress[0]++;
				
				if(progress[0]/size >= progress[1]) {
					long thisTime = System.currentTimeMillis();
					System.out.println(progress[0]+"/"+size+"\t.."+(thisTime-startTime)/1000+" sec");
					progress[1] += 0.1;
				}
				
				Record record = dupRecords.get(0);
				Trie frTrie = Trie.builder().addKeyword(record.frSequence).build();
				Trie rcTrie = Trie.builder().addKeyword(record.rcSequence).build();
				chr[0] = record.chr;
				
				SAMRecordIterator iterator = samReader.queryOverlapping(record.chr, record.start, record.end);
				int leastCount = 0;
	            while (iterator.hasNext()) {
	                SAMRecord samRecord = iterator.next();
	                // Process each SAM record
	                if(samRecord.getStart() <= record.start && samRecord.getEnd() >= record.end) {
		                String sequence = samRecord.getReadString();
		                Collection<Emit> emits = frTrie.parseText(sequence);
		                
		                if(emits.size() > 0) {
		                	record.frCnt++;
		                	leastCount++;
		                }
		                
		                emits = rcTrie.parseText(sequence);
		                if(emits.size() > 0) {
		                	record.rcCnt++;
		                	leastCount++;
		                }
	               }
	                
	            }
	            iterator.close();
	            
	            //
	            // in case of soft-clip, it can be zero because of unstable record range.
	            //////////////////////////////////// START SFOT-CLIP /////////////////////
	            if(leastCount == 0) {
	            	iterator = samReader.queryOverlapping(record.chr, record.start-100, record.end+100);
	            	while (iterator.hasNext()) {
		                SAMRecord samRecord = iterator.next();
		                // Process each SAM record
		                String sequence = samRecord.getReadString();
		                Collection<Emit> emits = frTrie.parseText(sequence);
		                
		                if(emits.size() > 0) {
		                	record.frCnt++;
		                	leastCount++;
		                }
		                
		                emits = rcTrie.parseText(sequence);
		                if(emits.size() > 0) {
		                	record.rcCnt++;
		                	leastCount++;
		                }
		                
		            }
		            iterator.close();
	            }
	            ////////////////////////////////////// END SFOT-CLIP /////////////////////
	            
	            for(int i=1; i<dupRecords.size(); i++) {
	            	dupRecords.get(i).frCnt = record.frCnt;
	            	dupRecords.get(i).rcCnt = record.rcCnt;
	            }
			});
			
			
			// for unmapped reads
			Trie trie = getTrie(task.records, Scan.mode);
			if(trie != null) {
				Hashtable<String, ArrayList<Record>> map = getMapBySequence(task.records);
				SAMRecordIterator iterator = samReader.queryUnmapped();
				while (iterator.hasNext()) {
	                SAMRecord samRecord = iterator.next();
	                // Process each SAM record
	                String sequence = samRecord.getReadString();
	                Collection<Emit> emits = trie.parseText(sequence);
					for(Emit emit : emits) {
						String matchedSequence = emit.getKeyword();
						ArrayList<Record> mRecords = map.get(matchedSequence);
						
						for(Record record : mRecords) {
							
							// only unmapped reads will be counted
							if(!record.chr.equalsIgnoreCase("*")) continue;
							
							if(record.frSequence.equalsIgnoreCase(matchedSequence)) {
								record.frCnt++;
							} else if(record.rcSequence.equalsIgnoreCase(matchedSequence)) {
								record.rcCnt++;
							} else {
								System.out.println("Out of scope!");
							}
						}
					}
	            }
	            iterator.close();
			}
            
		} catch(Exception e) {
			e.printStackTrace();
			System.out.println(chr[0]);
			System.exit(1);
		}
		long endTime = System.currentTimeMillis();
		System.out.println((endTime-startTime)/1000+" sec");
	}
	
}
