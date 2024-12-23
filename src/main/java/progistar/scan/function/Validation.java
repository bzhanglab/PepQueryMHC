package progistar.scan.function;

import progistar.scan.data.Constants;
import progistar.scan.data.Parameters;

public class Validation {

	private static String NUCLEOTIDE_REG_EXP = "A|C|T|G";
	private static String PEPTIDE_REG_EXP = "A|C|D|E|F|G|H|I|K|L|M|N|P|Q|R|S|T|V|W|Y";
	
	public static boolean checkValidSequence (String sequence) {
		boolean pass = true;
		
		if(Parameters.sequence.equalsIgnoreCase(Constants.SEQUENCE_NUCLEOTIDE)) {
			if(sequence.replaceAll(NUCLEOTIDE_REG_EXP, "").length() != 0) {
				pass = false;
			}
		} else if(Parameters.sequence.equalsIgnoreCase(Constants.SEQUENCE_PEPTIDE)) {
			if(sequence.replaceAll(PEPTIDE_REG_EXP, "").length() != 0) {
				pass = false;
			}
		}
		
		if(!pass) {
			System.out.println("Reject input: " + sequence);
		}
		
		return pass;
	}
}
