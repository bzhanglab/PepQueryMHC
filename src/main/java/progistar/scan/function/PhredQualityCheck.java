package progistar.scan.function;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.fastq.FastqRecord;
import progistar.scan.data.Parameters;
import progistar.scan.data.Phred;

public class PhredQualityCheck {

	
	public static boolean isPass (SAMRecord record, int start, int end) {
		String phredStr = record.getBaseQualityString().substring(start, end);
		return testByROI(phredStr);
	}
	
	public static boolean isPass (FastqRecord record, int start, int end) {
		String phredStr = record.getBaseQualityString().substring(start, end);
		return testByROI(phredStr);
	}
	
	/**
	 * quality check by single base.
	 * 
	 * @param phredStr
	 * @return
	 * @deprecated
	 */
	private static boolean testBySingleBase (String phredStr) {
		int len = phredStr.length();
		for(int i=0; i<len; i++) {
			int score = (int)phredStr.charAt(i);
			score -= 33; // Phred33
			
			if(score < Parameters.singleBaseThreshold) {
				return false;
			}
		}
		
		return true;
	}
	
	/**
	 * quality check by ROI unit. <br>
	 * discard if a probability of at least single base error is greater than ROIErrorThreshold. 
	 * 
	 * @param phredStr
	 * @return
	 */
	private static boolean testByROI (String phredStr) {
		if(Phred.getProbOfAtLeastOneError(phredStr) < Parameters.ROIErrorThreshold) {
			return true;
		}
		
		return false;
	}
}
