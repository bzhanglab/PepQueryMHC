package progistar.scan.data;

import htsjdk.samtools.util.Locatable;

public class BAMSRecord {

	public static String header = null;
	
	public String record;
	public String sequence;
	public String strand;
	public String location;
	
	public int start = 0;
	public int end = 0;
	public String chr = null;
	
	public int readCnt = 0;
	
	Locatable locatablePosition = null;
	
	public String getKey () {
		return (sequence+"_"+location+"_"+strand);
	}
}
