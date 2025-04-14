package progistar.scan.data;

import java.util.ArrayList;
import java.util.Collections;

public class Transcript {

	public int start;
	public int end;
	public String id;
	public ArrayList<String> tags;
	public boolean strand = true;
	public String warningTag = Constants.NULL;
	public ArrayList<Exon> exons = new ArrayList<Exon>();
	
	public void refine () {
		int size = exons.size();
		Collections.sort(exons);
		boolean hasCDS = false;
		for(int i=0; i<size; i++) {
			if(exons.get(i).feature == Constants.CDS) {
				hasCDS = true;
				break;
			}
		}
		if(hasCDS) {
			// refine Exon to CDS/UTR
			ArrayList<Exon> refines = new ArrayList<Exon>();
			
			// insert CDS only
			for(int i=0; i<size; i++) {
				if(exons.get(i).feature == Constants.CDS) {
					refines.add(exons.get(i));
				}
			}
			
			int cdsIdx = 0; 
			int cdsSize = refines.size();
			for(int i=0; i<size; i++) {
				Exon exon = exons.get(i);
				// pass CDS
				if(exon.feature == Constants.CDS) {
					continue;
				}
				
				Exon cds = refines.get(cdsIdx);
				// discard
				if(exon.start == cds.start && exon.end == cds.end) {
					cdsIdx++;
				}
				// upstream
				else if(exon.end < cds.start) {
					exon.feature = this.strand ? Constants.UTR5 : Constants.UTR3;
					refines.add(exon);
				} 
				// upstream overlap
				else if(exon.end == cds.end) {
					exon.feature = this.strand ? Constants.UTR5 : Constants.UTR3;
					exon.end = cds.start-1;
					cdsIdx++;
					refines.add(exon);
				}
				// downstream overlap
				else if(exon.start == cds.start) {
					exon.feature = this.strand ? Constants.UTR3 : Constants.UTR5;
					exon.start = cds.end+1;
					cdsIdx++;
					refines.add(exon);
				}
				// downstream
				else if(exon.start > cds.end) {
					exon.feature = this.strand ? Constants.UTR3 : Constants.UTR5;
					refines.add(exon);
				}
				// overlap
				else {
					Exon leftExon = new Exon();
					leftExon.start = exon.start;
					leftExon.end = cds.start-1;
					leftExon.feature = this.strand ? Constants.UTR5 : Constants.UTR3;
					
					Exon rightExon = new Exon();
					rightExon.start = cds.end+1;
					rightExon.end = exon.end;
					rightExon.feature = this.strand ? Constants.UTR3 : Constants.UTR5;
					
					refines.add(leftExon);
					refines.add(rightExon);
					cdsIdx++;
				}
				
				if(cdsIdx == cdsSize) {
					cdsIdx --;
				}
			}
			
			Collections.sort(refines);
			this.exons.clear();
			this.exons = refines;
		}
		
		size = exons.size();
		for(int i=0; i<size-1; i++) {
			Exon prevExon = exons.get(i);
			Exon nextExon = exons.get(i+1);
			
			if(nextExon.start - prevExon.end > 1) {
				Exon intron = new Exon();
				intron.start = prevExon.end+1;
				intron.end = nextExon.start-1;
				intron.feature = Constants.INTRON;
				exons.add(intron);
			}
			
		}
		
		Collections.sort(exons);
		for(int i=0; i<exons.size(); i++) {
			exons.get(i).idx = i;
		}
		
	}
	
	public Annotation getClassCode (SequenceRecord sRecord) {
		Annotation annotation	= new Annotation();
		annotation.transcript = this;
		
		String classCode = "Undefined";
		
		String[] locations	= sRecord.location.split("\\|");
		int[] qStarts	= new int[locations.length];
		int[] qEnds		= new int[locations.length];
		
		// is invalid location?
		for(int i=0; i<locations.length; i++) {
			// TODO: check validity
			String[] location = locations[i].split("\\:")[1].split("\\-");
			qStarts[i]	= Integer.parseInt(location[0]);
			qEnds[i]	= Integer.parseInt(location[1]);
		}
		
		// if the positions are out of bound to transcript ?
		boolean isIGR = false;
		for(int i=0; i<locations.length; i++) {
			if( (qStarts[i] < start || qStarts[i] > end) ||
				(qEnds[i] < start || qEnds[i] > end)) {
				isIGR = true;
				break;
			}
		}
		
		if(isIGR) {
			classCode = Constants.MARK_INTERGENIC;
			annotation.classCode = classCode;
			return annotation;
		}
		
		// NCDS, CDS, UTR5, UTR3, and Intron
		ArrayList<Exon> matchedExons = new ArrayList<Exon>();
		boolean isCDS		= true;
		boolean isUTR5		= false;
		boolean isUTR3		= false;
		boolean isNCDS		= false;
		boolean isIntron	= false;
		boolean isAS		= false;
		//System.out.println(sRecord.toString());
		for(Exon exon : exons) {
			// overlap
			for(int i=0; i<locations.length; i++) {
				int qStart	= qStarts[i];
				int qEnd	= qEnds[i];
				if( !((exon.start > qEnd) || (exon.end < qStart)) ) {
					matchedExons.add(exon);
					//System.out.println(this.id+": "+exon.feature+", "+exon.start+"-"+exon.end);
					if(exon.feature != Constants.CDS) {
						isCDS = false;
						
						switch(exon.feature) {
							case Constants.UTR5: isUTR5 = true; break;
							case Constants.UTR3: isUTR3 = true; break;
							case Constants.NCDS: isNCDS = true; break;
							case Constants.INTRON: isIntron = true; break;
							default: break;
						}
					}
				}
			}
		}
		
		// check AS
		if(locations.length != 1) {
			// if the number of reference exons is different from the size of locations
			// it implies that it should be mis-spliced.
			if(locations.length != matchedExons.size()) {
				isAS = true;
			} else {
				for(int i=0; i<locations.length; i++) {
					if( ((i%2 == 0) && (qEnds[i] == matchedExons.get(i).end)) || 
						((i%2 == 1) && (qStarts[i] == matchedExons.get(i).start))) {
						// match well
					} 
					else {
						isAS = true;
						break;
					}
				}
			}
		}
		
		// If CDS: In-frame? FS?
		// Else: ncRNA, UTR5, UTR3, Intron
		if(isCDS) {
			// find frame
			if(getFrameMark(qStarts[0], qEnds[qEnds.length-1]) == 0) {
				classCode = Constants.MARK_PC;
			} else {
				classCode = Constants.MARK_FS;
			}
			
		} else {
			// the worse is selected
			if(isIntron) {
				classCode = Constants.MARK_INTRON;
			} else if(isNCDS) {
				classCode = Constants.MARK_NCRNA;
			} else if(isUTR5) {
				classCode = Constants.MARK_UTR5;
			} else if(isUTR3) {
				classCode = Constants.MARK_UTR3;
			}
		}
		
		// if different strand?
		if(this.strand != (sRecord.strand.charAt(0) == '+') ) {
			classCode = Constants.MARK_ASRNA;
		}
		
		// AS check
		if(isAS) {
			classCode += ";" + Constants.MARK_AS;
		}
		
		annotation.classCode = classCode;
		
		return annotation;
	}
	
	public byte getFrameMark (int start, int end) {
		byte mark = Constants.FRAME_X;
		int pos = start;
		int cds = 0;
		int size = this.exons.size();
		
		if(this.strand) {
			pos = start;
			for(Exon exon : this.exons) {
				// inclusive
				if(exon.start <= pos && exon.end >= pos) {
					if(exon.feature == Constants.CDS) {
						cds += (pos - exon.start);
						mark = (byte) ( cds % 3);
					}
					break;
				}
				// exclusive
				else {
					if(exon.feature == Constants.CDS) {
						cds += (exon.end - exon.start + 1);
					}
				}
			}
		}
		
		else {
			pos = end;
			for(int i = size-1; i>=0; i--) {
				Exon exon = this.exons.get(i);
				// inclusive
				if(exon.start <= pos && exon.end >= pos) {
					if(exon.feature == Constants.CDS) {
						cds += (exon.end - pos);
						mark = (byte) ( cds % 3);
					}
					break;
				}
				// exclusive
				else {
					if(exon.feature == Constants.CDS) {
						cds += (exon.end - exon.start + 1);
					}
				}
			}
		}
		

		return mark;
	}
	
	/**
	 * Taking startLoci and endLoci of a region of mapped peptide. <br>
	 * return MARK_AS or MARK_CA. <br>
	 *
	 *
	 * @param startLoci
	 * @param endLoci
	 * @return
	 */
	public boolean isAS (ArrayList<Exon> matchedExons) {
		/**
		 * Think!
		 *
		 * EXON1 - INTRON1 - EXON2 - INTRON2 - EXON3
		 *
		 * 1) EXON1 - INTRON1 => is consecutive mapping?
		 * 2) EXON1 - EXON2 => is consecutive mapping?
		 * 3) EXON1 - EXON3 => okay. it must be AS
		 * 4) INTRON1 - INTRON2 => okay. it must be AS
		 *
		 */

		// non-consecutive mapping
		// this implies that
		// 1) EXON1 - INTRON1 => okay. it must be AS
		// 2) EXON1 - EXON2 => both ends are consecutive?
		// 3) EXON1 - EXON3 => okay. it must be AS
		// 4) INTRON1 - INTRON2 => okay. it must be AS
		// So! we just care about case 2. If not, it must be AS
		int size = matchedExons.size();
		// storing aBlock number for each start and end loci
		// if the locus resides on intergenic, store 0.

		boolean isAS = false;
		for(int i=0; i<size-1; i++) {
			Exon prevExon = matchedExons.get(i);
			Exon nextExon = matchedExons.get(i+1);
			
			if(prevExon.feature != Constants.INTRON && nextExon.feature != Constants.INTRON) {
				// is consecutive?
				if((prevExon.idx + 1) != nextExon.idx) {
					// junction variation?
					
					
				} else {
					isAS = true;
					break;
				}
			} else {
				isAS = true;
				break;
			}
		}
		

		return isAS;
	}
	
	
	public void printExons() {
		System.out.println(this.id+"\t"+this.strand);
		exons.forEach((e)->{
			String feature = ".";
			if(e.feature == Constants.CDS) {
				feature = "CDS";
			} else if(e.feature == Constants.NCDS) {
				feature = "NCDS";
			} else if (e.feature == Constants.UTR5) {
				feature = "5`-UTR";
			} else if (e.feature == Constants.UTR3) {
				feature = "3`-UTR";
			} else if(e.feature == Constants.INTRON) {
				feature = "IR";
			}
			
			System.out.println(feature+"("+e.start+"-"+e.end+")");
		});
	}
	
}
