package progistar.scan.data;

import htsjdk.samtools.util.Locatable;

public class Record {

	public static String header = null;
	
	public String record;
	public String frSequence;
	public String rcSequence;
	
	public int start = 0;
	public int end = 0;
	public String chr = null;
	public int genomicLociCount;
	
	public int frCnt = 0;
	public int rcCnt = 0;
	
	Locatable locatablePosition = null;
	
	public void reset() {
		this.frCnt = 0;
		this.rcCnt = 0;
	}
}
