package progistar.scan.data;

import java.util.ArrayList;
import java.util.Hashtable;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.ahocorasick.trie.Emit;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMTag;
import htsjdk.samtools.fastq.FastqRecord;
import progistar.scan.function.PhredQualityCheck;
import progistar.scan.function.Translator;
import progistar.scan.function.Utils;

public class LocationInformation {
	public static final Pattern EACH_CIGAR_REGEX = Pattern.compile("([0-9]+)([MINDSHPX=])");
	public static final Pattern EACH_MD_REGEX = Pattern.compile("(([0-9]+)|([A-Z]+|\\^[A-Z]+))");
	
	public String inputSequence;
	
	public String location = Constants.NULL;
	public String mutation = Constants.NULL;
	public String obsNucleotide = Constants.NULL;
	public String refNucleotide = Constants.NULL;
	public String obsPeptide = Constants.NULL;
	public String refPeptide = Constants.NULL;
	public Hashtable<String, Long> readCounts = new Hashtable<String, Long>();
	public char strand = Constants.NULL.charAt(0);
	
	public String getKey () {
		return location+"\t"+strand+"\t"+obsNucleotide+"\t"+refNucleotide;
	}
	
	public String getRes () {
		if(Parameters.isSingleCellMode) {
			StringBuilder str = new StringBuilder(location+"\t"+mutation+"\t"+strand+"\t"+obsPeptide+"\t"+obsNucleotide+"\t"+refNucleotide);
			// write raw read counts
			for(String barcodeId : BarcodeTable.barcodeIds) {
				Long read = readCounts.get(barcodeId);
				if(read == null) {
					read = 0L;
				}
				str.append("\t").append(read);
			}
			// write RPHTs
			for(String barcodeId : BarcodeTable.barcodeIds) {
				Long read = readCounts.get(barcodeId);
				if(read == null) {
					read = 0L;
				}
				str.append("\t").append(Utils.getRPHT(read, barcodeId));
			}
			return str.toString();
		} else {
			Long read = readCounts.get(Constants.DEFAULT_BARCODE_ID);
			return location+"\t"+mutation+"\t"+strand+"\t"+obsPeptide+"\t"+obsNucleotide+"\t"+refNucleotide+"\t"+read+"\t"+Utils.getRPHM((double)read, Constants.DEFAULT_BARCODE_ID);
		}
	}
	
	public long getTotalReads () {
		long sum = 0;
		if(Parameters.isSingleCellMode) {
			for(String barcodeId : BarcodeTable.barcodeIds) {
				Long read = readCounts.get(barcodeId);
				if(read == null) {
					read = 0L;
				}
				sum += read;
			}
		} else {
			Long read = readCounts.get(Constants.DEFAULT_BARCODE_ID);
			sum += read;
		}
		
		return sum;
	}
	
	public void calMetaInfo () {
		this.calMutation();
		assert this.obsNucleotide != null;
		
		if(strand == '+') {
			this.obsPeptide = Translator.translation(this.obsNucleotide, 0);
		} else if(strand == '-') {
			this.obsPeptide = Translator.translation(Translator.getReverseComplement(this.obsNucleotide), 0);
		}
		// if a strand is not defined, it should save forward strand nucleotide
		else if(strand == Constants.NULL.charAt(0)) {
			this.obsPeptide = Translator.translation(this.obsNucleotide, 0);
		}
	}
	
	public void calMutation () {
		String[] locations = location.split("\\|");
		int mPos =0;
		ArrayList<Mutation> mutations = new ArrayList<Mutation>();
		for(int i=0; i<locations.length; i++) {
			String location = locations[i];
			
			
			try {
				// unmapped location
				if(location.equalsIgnoreCase(Constants.NULL)) {
					this.mutation = Constants.NULL;
				} else {
					String chr = location.split("\\:")[0];
					String[] fields = location.split("\\:")[1].split("\\-");
					int start = Integer.parseInt(fields[0]);
					int end = Integer.parseInt(fields[1]);
					
					for(int j=start; j<=end; j++) {
						// mutation
						if(this.obsNucleotide.charAt(mPos) != this.refNucleotide.charAt(mPos)) {
							Mutation mutation = new Mutation();
							mutation.altSeq = this.obsNucleotide.charAt(mPos)+"";
							mutation.chrName = chr;
							mutation.refSeq = this.refNucleotide.charAt(mPos)+"";
							mutation.genomicPosition = j;
							// softclip
							if(this.refNucleotide.charAt(mPos) == '*') {
								mutation.type = Constants.CLP;
							} 
							// insertion
							else if(this.refNucleotide.charAt(mPos) == '.') {
								j--;
								mutation.genomicPosition = j;
								mutation.type = Constants.INS;
							}
							// deletion
							else if(this.obsNucleotide.charAt(mPos) == '.') {
								mutation.type = Constants.DEL;
							}
							// snp
							else {
								mutation.type = Constants.SNP;
							}
							mutations.add(mutation);
						}
						mPos++;
					}
				}
			}catch(Exception e) {
				e.printStackTrace();
				
				System.out.println(location);
			}
			
		}
		
		if(mutations.size() == 0) {
			this.mutation = Constants.NULL;
		} else {
			this.mutation = "";
			for(int i=0; i<mutations.size()-1; i++) {
				Mutation prevM = mutations.get(i);
				Mutation nextM = mutations.get(i+1);
				
				if(prevM.type == Constants.SNP) {
					continue;
				}
				
				if(prevM.type == nextM.type) {
					prevM.refSeq += nextM.refSeq;
					prevM.altSeq += nextM.altSeq;
					mutations.remove(i+1);
					i--;
				}
			}
			
			for(int i=0; i<mutations.size(); i++) {
				// lower to upper
				mutations.get(i).refSeq = mutations.get(i).refSeq.toUpperCase();
				if(mutation.length() != 0) {
					mutation += "|";
				}
				mutation += mutations.get(i).toString();
					
			}
		}
	}
	

	/**
	 * 
	 * Return genomic location of a given sequence.<br>
	 * If there is no available position in the SAMRecord, return ".". <br>
	 * If the sequence contains low quality position-specific Phred score, than it return null.
	 * 
	 * @param samRecord
	 * @param emit
	 * @param frame
	 * @return
	 */
	public static LocationInformation getMatchedLocation (SAMRecord samRecord, Emit emit, int frame, char strand) {
		boolean isMD = true;
		Object mdTag = samRecord.getAttribute(SAMTag.MD);
		String barcodeId = BarcodeTable.getBarcodeFromBam(samRecord);
		
		if(mdTag == null) {
			isMD = false;
		}
		
		LocationInformation lInfo = new LocationInformation();
		lInfo.strand = strand;
		lInfo.readCounts.put(barcodeId, 1L);
		
		String nucleotide = samRecord.getReadString();
		String reference = getRefSequence(samRecord);
		
		StringBuilder obsSequenceGivenRegion = new StringBuilder();
		StringBuilder refSequenceGivenRegion = new StringBuilder();
		
		int startPos = emit.getStart();
		int endPos = emit.getEnd()+1;
		
		if(Parameters.sequence.equalsIgnoreCase(Constants.SEQUENCE_PEPTIDE)) {
			startPos = (startPos) * 3 + frame;
			endPos = (endPos) * 3 + frame;
		}
		
		if(strand == '-') {
			int len = nucleotide.length();
			int tmp = len - startPos;
			startPos = len - endPos;
			endPos = tmp;
		}
		
		// [startPos, endPos) zero-based
		
		// check quality
		boolean isPass = PhredQualityCheck.isPass(samRecord, startPos, endPos);
		if(!isPass) {
			return null;
		}
		
		String cigarStr = samRecord.getCigarString();
		// unmapped reads
		if(cigarStr.equalsIgnoreCase("*")) {
			// 
			lInfo.location = Constants.NULL;
			lInfo.strand = Constants.NULL.charAt(0);
			lInfo.obsNucleotide = nucleotide.substring(startPos, endPos);
			lInfo.refNucleotide = Constants.NULL;//reference.substring(startPos, endPos);
			
			if(strand == '-') {
				lInfo.obsNucleotide = Translator.getReverseComplement(lInfo.obsNucleotide);
				lInfo.refNucleotide = Constants.NULL;//Translator.getReverseComplement(lInfo.refNucleotide);
			}
			
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
		
		/**
		 * Out of range (e.g., ChrM:-3-15)<br>
		 * Softclip at the end of sequence and positioned at the near start location in a chromosome.
		 * 
		 * 
		 */
		boolean isOutOfRange = false;
		for(int i=0; i<startLocations.size(); i++) {
			int start = startLocations.get(i);
			int end = endLocations.get(i);
			
			if(start <= 0 || end <= 0) {
				isOutOfRange = true;
			}
			
			if(locations.length() != 0) {
				locations.append("|");
			}
			locations.append(chr+":"+start+"-"+end);
		}
		
		if(isOutOfRange) {
			locations.setLength(0);
			locations.append(Constants.NULL);
		}
		
		lInfo.location = locations.toString();
		lInfo.obsNucleotide = obsSequenceGivenRegion.toString();
		lInfo.refNucleotide = refSequenceGivenRegion.toString();
		return lInfo;
	}
	
	public static LocationInformation getMatchedLocation (FastqRecord fastqRecord, Emit emit, int frame, char strand) {
		String barcodeId = BarcodeTable.getBarcodeFromFASTQ(fastqRecord);
		
		LocationInformation lInfo = new LocationInformation();
		lInfo.strand = strand;
		lInfo.readCounts.put(barcodeId, 1L);
		
		String nucleotide = fastqRecord.getReadString();
		
		int startPos = emit.getStart();
		int endPos = emit.getEnd()+1;
		
		if(Parameters.sequence.equalsIgnoreCase(Constants.SEQUENCE_PEPTIDE)) {
			startPos = (startPos) * 3 + frame;
			endPos = (endPos) * 3 + frame;
		}
		
		if(strand == '-') {
			int len = nucleotide.length();
			int tmp = len - startPos;
			startPos = len - endPos;
			endPos = tmp;
		}
		
		// [startPos, endPos) zero-based
		
		// check quality
		boolean isPass = PhredQualityCheck.isPass(fastqRecord, startPos, endPos);
		if(!isPass) {
			return null;
		}
		
		// 
		lInfo.location = Constants.NULL;
		lInfo.strand = Constants.NULL.charAt(0);
		lInfo.obsNucleotide = nucleotide.substring(startPos, endPos);
		lInfo.refNucleotide = Constants.NULL;
		
		if(strand == '-') {
			lInfo.obsNucleotide = Translator.getReverseComplement(lInfo.obsNucleotide);
			lInfo.refNucleotide = Constants.NULL;//Translator.getReverseComplement(lInfo.refNucleotide);
		}
		
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
