package progistar.scan.function;

import java.util.regex.Pattern;

import progistar.scan.data.Constants;
import progistar.scan.data.Parameters;

public class Validation {

	private static final Pattern NUCLEOTIDE_REG_EXP = Pattern.compile("A|C|T|G");
	private static final Pattern PEPTIDE_REG_EXP = Pattern.compile("A|C|D|E|F|G|H|I|K|L|M|N|P|Q|R|S|T|V|W|Y");

	public static boolean checkValidSequence (String sequence) {
		boolean pass = true;

		if(Parameters.sequence.equalsIgnoreCase(Constants.SEQUENCE_NUCLEOTIDE)) {
			if(NUCLEOTIDE_REG_EXP.matcher(sequence).replaceAll("").length() != 0) {
				pass = false;
			}
		} else if(Parameters.sequence.equalsIgnoreCase(Constants.SEQUENCE_PEPTIDE)) {
			if(PEPTIDE_REG_EXP.matcher(sequence).replaceAll("").length() != 0) {
				pass = false;
			}
		}
		
		if(!pass) {
			System.out.println("Reject input: " + sequence);
		}
		
		return pass;
	}
}
