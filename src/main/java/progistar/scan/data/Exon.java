package progistar.scan.data;

/**
 * Annotation Block <br>
 * 
 * @author gistar
 *
 */
public class Exon implements Comparable<Exon>{

	public int idx;
	public int start;
	public int end;
	// features:
	// CDS, UTR5, UTR3, NCDS, INTRON, INTERGENIC in Constants class
	public byte feature;
	
	@Override
	public int compareTo(Exon o) {
		if(this.start < o.start) {
			return -1;
		}else if(this.start > o.start) {
			return 1;
		}
		
		return 0;
	}
	
	public int getLength () {
		return (end - start + 1);
	}
}
