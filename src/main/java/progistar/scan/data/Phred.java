package progistar.scan.data;

public class Phred {

	private static double[] ERROR_PROB_TABLE = new double[101];
	private static double[] CORRECT_PROB_TABLE = new double[101];
	private static double[] CORRECT_PHRED_TABLE = new double[101];
	
	private static double convertToProb (double phred) {
		return Math.pow(10, -phred/10);
	}
	
	private static double convertToPhred (double prob) {
		return -10 * Math.log10(prob);
	}
	
	public static void loadTable () {
		for(int phred=0; phred<ERROR_PROB_TABLE.length; phred++) {
			ERROR_PROB_TABLE[phred] = convertToProb(phred);
			CORRECT_PROB_TABLE[phred] = 1 - ERROR_PROB_TABLE[phred];
			CORRECT_PHRED_TABLE[phred] = phred == 0 ? 0 : convertToPhred(CORRECT_PROB_TABLE[phred]);
		}
	}
	
	/**
	 * Use ROI
	 * @param phredStr
	 * @return
	 */
	private static double getProbOfAllCorrected (String phredStr) {
		double sumOfCorrectPhreds = 0;
		int len = phredStr.length();
		for(int i=0; i<len; i++) {
			int phred = (int)phredStr.charAt(i);
			phred -= 33; // Phred33
			
			sumOfCorrectPhreds += CORRECT_PHRED_TABLE[phred];
		}
		
		return convertToProb(sumOfCorrectPhreds);
	}
	
	public static double getProbOfAtLeastOneError (String phredStr) {
		return 1 - getProbOfAllCorrected(phredStr);
	}
}
