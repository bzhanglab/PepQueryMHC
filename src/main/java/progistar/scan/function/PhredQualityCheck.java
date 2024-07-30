package progistar.scan.function;

import htsjdk.samtools.SAMRecord;
import progistar.scan.run.Scan;

public class PhredQualityCheck {

	
	public static boolean isPass (SAMRecord record, int start, int end) {
		boolean isPass = true;
		String phredQuality = record.getBaseQualityString();
		
		for(int i=start; i<end; i++) {
			int score = (int)phredQuality.charAt(i);
			score -= 33; // Phred33
			
			if(score < Scan.phredThreshold) {
				isPass = false;
				break;
			}
		}
		
		return isPass;
	}
}
