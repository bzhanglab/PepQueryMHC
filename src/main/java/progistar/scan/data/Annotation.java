package progistar.scan.data;

import java.util.Hashtable;
import java.util.LinkedList;

public class Annotation implements Comparable<Annotation> {

	public String classCode;
	public double totalPenalty;
	public Transcript transcript;
	/**
	 * Warning tag is changed under the conditions:
	 * 1) If there is no transcript information:
	 * 	a) No location information provided.
	 * 	b) No matched reference name provided.
	 * 
	 * 2) If there is a matched transcript information:
	 * 	a) borrow a warning tag from transcript class.
	 */
	public String warningTag = Constants.NULL;
	public Gene gene;
	
	public static LinkedList<Annotation> removeRedundancy (LinkedList<Annotation> annotations) {
		Hashtable<String, Annotation> saveUnique = new Hashtable<String, Annotation>();
		LinkedList<Annotation> newAnnotations = new LinkedList<Annotation> ();
		annotations.forEach((a)->{
			saveUnique.put(a.key(), a);
		});
		
		saveUnique.forEach((k, a)->{
			newAnnotations.add(a);
		});
		
		return newAnnotations;
	}
	
	@Override
	/**
	 * Lower is better
	 */
	public int compareTo(Annotation o) {
		if(this.totalPenalty < o.totalPenalty) {
			return -1;
		} else if(this.totalPenalty > o.totalPenalty) {
			return 1;
		}
		return 0;
	}
	
	public String key () {
		String geneId = getGeneId();
		String geneName = getGeneName();
		String geneType = getGeneType();
		String strand = getStrand();
		
		return geneId + "_" + geneName + "_" + geneType + "_" + strand + "_" + classCode;
	}
	
	public String toString () {
		String geneId = getGeneId();
		String geneName = getGeneName();
		String geneType = getGeneType();
		String strand = getStrand();
		String tId = getTranscriptId();
		
		return geneId + "(" + geneName + ")" + ", " + tId + ", "+ geneType +", "+ strand +", "+classCode;
	}
	
	public String getGeneName () {
		return gene == null ? Constants.NULL : gene.name;
	}
	
	public String getGeneId () {
		return gene == null ? Constants.NULL : gene.id;
	}
	
	public String getGeneType () {
		return gene == null ? Constants.NULL : gene.type;
	}
	
	public String getStrand () {
		return transcript == null ? Constants.NULL : (transcript.strand ? "+" : "-");
	}
	
	public String getClassCode () {
		return classCode;
	}
	
	public String getTranscriptId () {
		return transcript == null ? Constants.NULL : transcript.id;
	}
	
	public String getWarningTag () {
		return this.warningTag;
	}
	
	public void complete () {
		this.calPenalty();
		// nothing to do if this is unknown.
		if(this.classCode.contains(Parameters.MARK_UNKNOWN)) {
			return;
		}
		
		
		// region
		/// the input class code should be naive class code.
		String representative = "UNDEF";
		String structuralAnnotations = "";
		int currentPriority = Integer.MAX_VALUE;
		// Class I: Exonic regions such as FS
		if(this.classCode.contains(Parameters.MARK_IF)) {
			if(Constants.ANNOTATION_PRIORITY_IF < currentPriority) {
				representative = Parameters.MARK_IF;
				currentPriority = Constants.ANNOTATION_PRIORITY_IF;
			}
		}
		if(this.classCode.contains(Parameters.MARK_UTR3)) {
			if(Constants.ANNOTATION_PRIORITY_3UTR < currentPriority) {
				representative = Parameters.MARK_UTR3;
				currentPriority = Constants.ANNOTATION_PRIORITY_3UTR;
			}
		}
		if(this.classCode.contains(Parameters.MARK_UTR5)) {
			if(Constants.ANNOTATION_PRIORITY_5UTR < currentPriority) {
				representative = Parameters.MARK_UTR5;
				currentPriority = Constants.ANNOTATION_PRIORITY_5UTR;
			}
		}
		if(this.classCode.contains(Parameters.MARK_OOF)) {
			if(Constants.ANNOTATION_PRIORITY_OOF < currentPriority) {
				representative = Parameters.MARK_OOF;
				currentPriority = Constants.ANNOTATION_PRIORITY_OOF;
			}
		}
		if(this.classCode.contains(Parameters.MARK_NCRNA)) {
			if(Constants.ANNOTATION_PRIORITY_NCRNA < currentPriority) {
				representative = Parameters.MARK_NCRNA;
				currentPriority = Constants.ANNOTATION_PRIORITY_NCRNA;
			}
		}
		
		// Class II: IR
		if(this.classCode.contains(Parameters.MARK_INTRON)) {
			if(Constants.ANNOTATION_PRIORITY_IR < currentPriority) {
				representative = Parameters.MARK_INTRON;
				currentPriority = Constants.ANNOTATION_PRIORITY_IR;
			}
		}
		// Class III: asRNA
		if(this.classCode.contains(Parameters.MARK_ASRNA)) {
			if(Constants.ANNOTATION_PRIORITY_ASRNA < currentPriority) {
				representative = Parameters.MARK_ASRNA;
				currentPriority = Constants.ANNOTATION_PRIORITY_ASRNA;
			}
		}
		// Class IV: IGR
		if(this.classCode.contains(Parameters.MARK_INTERGENIC)) {
			if(Constants.ANNOTATION_PRIORITY_IGR < currentPriority) {
				representative = Parameters.MARK_INTERGENIC;
				currentPriority = Constants.ANNOTATION_PRIORITY_IGR;
			}
		}
		
		// Class V: Structural annotation
		if(this.classCode.contains(Parameters.MARK_ES)) {
			if(structuralAnnotations.length() > 0) {
				structuralAnnotations += ";";
			}
			structuralAnnotations += Parameters.MARK_ES;
		}
		if(this.classCode.contains(Parameters.MARK_EE)) {
			if(structuralAnnotations.length() > 0) {
				structuralAnnotations += ";";
			}
			structuralAnnotations += Parameters.MARK_EE;
		}
		
		// add region + struct
		if(structuralAnnotations.length() > 0) {
			representative += ";" + structuralAnnotations;
		}
		// System.out.println(this.classCode +" => "+representative);
		this.classCode = representative;
	}
	
	private void calPenalty () {
		
		// Unknown Class
		if(this.classCode.contains(Parameters.MARK_UNKNOWN)) {
			this.totalPenalty = Constants.PENALTY_UNMAP;
			return;
		}
		
		double pClassI = 0;
		double pClassII = 0;
		double pClassIII = 0;
		double pClassIV = 0;
		double pClassV = 0;
		
		// Class I: Exonic regions such as FS
		if(this.classCode.contains(Parameters.MARK_UTR3)) {
			pClassI = Math.max(pClassI, Constants.PENALTY_3UTR);
		}
		if(this.classCode.contains(Parameters.MARK_UTR5)) {
			pClassI = Math.max(pClassI, Constants.PENALTY_5UTR);
		}
		if(this.classCode.contains(Parameters.MARK_OOF)) {
			pClassI = Math.max(pClassI, Constants.PENALTY_OOF);
		}
		if(this.classCode.contains(Parameters.MARK_NCRNA)) {
			pClassI = Math.max(pClassI, Constants.PENALTY_NCRNA);
		}
		
		// Class II: IR
		if(this.classCode.contains(Parameters.MARK_INTRON)) {
			pClassII = Constants.PENALTY_IR;
		}
		// Class III: asRNA
		if(this.classCode.contains(Parameters.MARK_ASRNA)) {
			pClassIV = Constants.PENALTY_ASRNA;
		}
		// Class IV: IGR
		if(this.classCode.contains(Parameters.MARK_INTERGENIC)) {
			pClassIII = Constants.PENALTY_IGR;
		}
		
		
		// Class V: SV
		if(this.classCode.contains(Parameters.MARK_ES)) {
			pClassV += Constants.PENALTY_ES;
		}
		if(this.classCode.contains(Parameters.MARK_EE)) {
			pClassV += Constants.PENALTY_EE;
		}
		
		
		this.totalPenalty = pClassI + pClassII + pClassIII + pClassIV + pClassV;
		// check warning code
		if(!this.getWarningTag().equalsIgnoreCase(Constants.NULL)) {
			this.totalPenalty = Constants.PENALTY_WARNING;
			
			if(Parameters.verbose) {
				System.out.println("Detect warning tags: "+this.transcript.warningTag +" in "+this.transcript.id);
				System.out.println("This annotation will get the maximum penalty.");
			}
		}
	}
	
}
