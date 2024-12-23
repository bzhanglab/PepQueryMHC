package progistar.scan.data;

import java.util.Hashtable;
import java.util.LinkedList;

public class Annotation implements Comparable<Annotation> {

	public String classCode;
	public int penalty;
	public Transcript transcript;
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
		if(this.penalty < o.penalty) {
			return -1;
		} else if(this.penalty > o.penalty) {
			return 1;
		}
		return 0;
	}
	
	public String key () {
		String geneId = gene == null ? Constants.NULL : gene.id;
		String geneName = gene == null ? Constants.NULL : gene.name;
		String geneType = gene == null ? Constants.NULL : gene.type;
		String strand = transcript == null ? Constants.NULL : (transcript.strand ? "+" : "-");
		
		return geneId + "@" + geneName + "@" + geneType + "@" + strand + "@" + classCode;
	}
	
	public String toString () {
		String geneId = gene == null ? Constants.NULL : gene.id;
		String geneName = gene == null ? Constants.NULL : gene.name;
		String geneType = gene == null ? Constants.NULL : gene.type;
		String strand = transcript == null ? Constants.NULL : (transcript.strand ? "+" : "-");
		String tId = transcript == null ? Constants.NULL : transcript.id;
		
		
		return geneId + "(" + geneName + ")" + ", " + tId + ", "+ geneType +", "+ strand +", "+classCode;
	}
	
	public void calPenalty () {
		this.penalty = 0;
		if(this.classCode.contains(Constants.MARK_AS)) {
			this.penalty += Constants.PENALTY_AS;
		}
		if(this.classCode.contains(Constants.MARK_NCRNA)) {
			this.penalty += Constants.PENALTY_ncRNA;
		}
		if(this.classCode.contains(Constants.MARK_ASRNA)) {
			this.penalty += Constants.PENALTY_asRNA;
		}
		if(this.classCode.contains(Constants.MARK_INTRON)) {
			this.penalty += Constants.PENALTY_IR;
		}
		if(this.classCode.contains(Constants.MARK_INTERGENIC)) {
			this.penalty += Constants.PENALTY_IGR;
		}
		if(this.classCode.contains(Constants.MARK_UTR3)) {
			this.penalty += Constants.PENALTY_3UTR;
		}
		if(this.classCode.contains(Constants.MARK_UTR5)) {
			this.penalty += Constants.PENALTY_5UTR;
		}
		if(this.classCode.contains(Constants.MARK_UNKNOWN)) {
			this.penalty += Constants.PENALTY_UNMAP;
		}
	}
	
}
