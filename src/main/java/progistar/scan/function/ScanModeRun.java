package progistar.scan.function;

import java.io.File;
import java.util.ArrayList;
import java.util.Collection;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.ahocorasick.trie.Emit;
import org.ahocorasick.trie.Trie;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SAMTag;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import progistar.scan.data.SequenceRecord;
import progistar.scan.data.Constants;
import progistar.scan.data.LocationInformation;
import progistar.scan.run.Scan;
import progistar.scan.run.Task;

public class ScanModeRun {
	public static final Pattern EACH_CIGAR_REGEX = Pattern.compile("([0-9]+)([MINDSHPX=])");
	public static final Pattern EACH_MD_REGEX = Pattern.compile("(([0-9]+)|([A-Z]+|\\^[A-Z]+))");
	
	
	public static void runScanMode (Task task) {
		if(task.type == Constants.TYPE_SCAN_MODE_TASK) {
			System.out.println(task.chrName+":"+task.start+"-"+task.end);
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
		System.out.println("Task"+task.taskIdx+" "+(endTime-startTime)/1000+" sec");
	}
	
	
	private static void find (SAMRecordIterator iterator, Trie trie, Task task) {
		int count = 0;
		while (iterator.hasNext()) {
            SAMRecord samRecord = iterator.next();
            count ++;
            
            boolean isPass = false;
            // if the task is for mapped reads
            // only reads with below that genomic start are retrieved
            if(task.readType == Constants.MAPPED_READS) {
            	if( !(samRecord.getAlignmentStart() >= task.start && 
            			samRecord.getAlignmentStart() < task.end) ) {
            		isPass = true;
            	}
            	
            	if(Scan.count.equalsIgnoreCase(Constants.COUNT_PRIMARY) && samRecord.isSecondaryAlignment()) {
            		isPass = true;
            	}
            }
            // if the task is for unmapped reads
            // only process reads given range
            else {
            	if(count < task.start || count > task.end) {
            		isPass = true;
            	}
            }

            
            if(isPass) {
            	continue;
            }
            
            // increase processed reads if and only if
            // a read is primary and mapped to a genomic region.
            if(!samRecord.isSecondaryAlignment() && task.readType == Constants.MAPPED_READS) {
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
        				LocationInformation matchedLocation = getMatchedLocation(samRecord, emit, 0, strandChar);
        				matchedLocation.inputSequence = emit.getKeyword();
        				
        				if(task.locTable.putLocation(matchedLocation)) {
        					matchedLocation.calMetaInfo();
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
            				LocationInformation matchedLocation = getMatchedLocation(samRecord, emit, fr, strandChar);
            				matchedLocation.inputSequence = emit.getKeyword();
            				
            				if(task.locTable.putLocation(matchedLocation)) {
            					matchedLocation.calMetaInfo();
            				}
            			}
            		}
            	}
            }
            
        }
        iterator.close();
	}
	
	
	/**
	 * 
	 * Return genomic location of a given sequence.<br>
	 * If there is no available position in the SAMRecord, return ".".
	 * 
	 * @param samRecord
	 * @param emit
	 * @param frame
	 * @return
	 */
	private static LocationInformation getMatchedLocation (SAMRecord samRecord, Emit emit, int frame, char strand) {
		boolean isMD = true;
		Object mdTag = samRecord.getAttribute(SAMTag.MD);
		
		if(mdTag == null) {
			isMD = false;
		}
		
		LocationInformation lInfo = new LocationInformation();
		lInfo.strand = strand;
		
		String nucleotide = samRecord.getReadString();
		String reference = getRefSequence(samRecord);
		
		StringBuilder obsSequenceGivenRegion = new StringBuilder();
		StringBuilder refSequenceGivenRegion = new StringBuilder();
		
		// [startPos, endPos] zero-based
		int startPos = emit.getStart();
		int endPos = emit.getEnd()+1;
		
		if(Scan.sequence.equalsIgnoreCase(Constants.SEQUENCE_PEPTIDE)) {
			startPos = (startPos) * 3 + frame;
			endPos = (endPos) * 3 + frame;
		}
		
		if(strand == '-') {
			int len = nucleotide.length();
			int tmp = len - startPos;
			startPos = len - endPos;
			endPos = tmp;
		}
		
		String cigarStr = samRecord.getCigarString();
		// unmapped reads
		if(cigarStr.equalsIgnoreCase("*")) {
			lInfo.location = Constants.NULL;
			lInfo.obsNucleotide = nucleotide.substring(startPos, endPos);
			lInfo.refNucleotide = reference.substring(startPos, endPos);
			return lInfo;
		}
		
		String chr = samRecord.getReferenceName();
		
		
		Matcher matcher = EACH_CIGAR_REGEX.matcher(cigarStr);
		ArrayList<Integer> startLocations = new ArrayList<Integer>();
		ArrayList<Integer> endLocations = new ArrayList<Integer>();
		int baseStartPos = samRecord.getUnclippedStart(); // one-based
		
		int startGenomicPosition = -1;
		int endGenomicPosition = -1;
		
		int seqPos = -1;
		int refPos = -1;
		int gPos = baseStartPos - 1;
		
		while (matcher.find()) {
			int markerSize = Integer.parseInt(matcher.group(1));
			char operation = matcher.group(2).charAt(0);
			
			if(operation == 'N') {
				if(startGenomicPosition != -1 && endGenomicPosition != -1) {
					startLocations.add(startGenomicPosition);
					endLocations.add(endGenomicPosition);
				}
				
				startGenomicPosition = -1;
				endGenomicPosition = -1;
				gPos += markerSize;
			} else if(operation == 'M' || operation == 'D' || operation == 'S' || operation == 'I') {
				for(int i=0; i<markerSize; i++) {
					if(operation == 'D') {
						// do not consume sequence position
						if(isMD) {
							refPos++;
						}
						gPos++;
					} else if (operation == 'I') {
						seqPos ++;
						refPos ++;
					} else {
						gPos++;
						seqPos ++;
						refPos ++;
					}
					
					try {
						if(seqPos >= startPos && seqPos < endPos) {
							if(startGenomicPosition == -1 && operation != 'I') {
								startGenomicPosition = gPos;
							}
							endGenomicPosition = gPos;
							
							if(operation == 'D' ) {
								obsSequenceGivenRegion.append(".");
								if(isMD) {
									refSequenceGivenRegion.append(reference.charAt(refPos));
								} else {
									refSequenceGivenRegion.append(".");
								}
							} else if (operation == 'I') {
								obsSequenceGivenRegion.append(nucleotide.charAt(seqPos));
								refSequenceGivenRegion.append(reference.charAt(refPos));
							} else {
								obsSequenceGivenRegion.append(nucleotide.charAt(seqPos));
								refSequenceGivenRegion.append(reference.charAt(refPos));
							}
						}
					}catch(Exception e) {
						e.printStackTrace();
						System.out.println(emit.getKeyword());
						System.out.println(samRecord.getSAMString());
						
						System.exit(-1);
					}
					
				}
			} else {
				System.out.println(operation);
			}
	    }
		
		if(startGenomicPosition != -1 && endGenomicPosition != -1) {
			startLocations.add(startGenomicPosition);
			endLocations.add(endGenomicPosition);
		}
		
		StringBuilder locations = new StringBuilder();
		for(int i=0; i<startLocations.size(); i++) {
			int start = startLocations.get(i);
			int end = endLocations.get(i);
			if(locations.length() != 0) {
				locations.append("|");
			}
			locations.append(chr+":"+start+"-"+end);
		}
		
		lInfo.location = locations.toString();
		lInfo.obsNucleotide = obsSequenceGivenRegion.toString();
		lInfo.refNucleotide = refSequenceGivenRegion.toString();
		return lInfo;
	}
	
	private static String getRefSequence (SAMRecord samRecord) {
		StringBuilder refSequence = new StringBuilder();
		String obsSequence = samRecord.getReadString();
		refSequence.append(obsSequence);
		Object mdTag = samRecord.getAttribute(SAMTag.MD);
		
		if(mdTag == null) {
			return refSequence.toString();
		}
		
		// find softclip
		Matcher matcher = EACH_CIGAR_REGEX.matcher(samRecord.getCigarString());
		int pos = -1;
		while (matcher.find()) {
			int markerSize = Integer.parseInt(matcher.group(1));
			char operation = matcher.group(2).charAt(0);
			
			if(operation == 'S') {
				for(int i=0; i<markerSize; i++) {
					pos++;
					refSequence.setCharAt(pos, '*');
				}
			} else if(operation == 'I') {
				for(int i=0; i<markerSize; i++) {
					pos++;
					refSequence.setCharAt(pos, '.');
				}
			} else if(operation == 'M') {
				pos += markerSize;
			}
		}

		// MD parsing
		Matcher mdMatcher = EACH_MD_REGEX.matcher((String) mdTag);
		pos = 0;
		while(mdMatcher.find()) {
			String md = mdMatcher.group();
			char sign = md.charAt(0);
			
			// match size
			if(Character.isDigit(sign)) {
				pos += Integer.parseInt(md);
			} 
			// nt change
			else if(Character.isAlphabetic(sign)) {
				for(int i=0; i<md.length(); i++) {
					pos++;
					
					int mPos = 0;
					for(int j=0; j<refSequence.length(); j++) {
						if(refSequence.charAt(j) != '*' && refSequence.charAt(j) != '.') {
							mPos++;
							
							if(mPos == pos) {
								refSequence.setCharAt(j, Character.toLowerCase(md.charAt(i)));
							}
						}
					}
				}
			} 
			// deletion sequence
			else if(sign == '^') {
				int mPos = 0;
				for(int j=0; j<refSequence.length(); j++) {
					if(refSequence.charAt(j) != '*' && refSequence.charAt(j) != '.') {
						mPos++;
						
						if(mPos == pos) {
							refSequence.insert(j+1, md.substring(1).toLowerCase());
							break;
						}
					}
				}
			}
		}
		
		
		return refSequence.toString();
	}
	

}











