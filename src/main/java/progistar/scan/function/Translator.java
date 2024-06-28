package progistar.scan.function;

import progistar.scan.data.Codon;
import progistar.scan.data.Constants;

public class Translator {
	
	public static String translation (String nucleotides, int frame) {
		StringBuilder peptides = new StringBuilder();
		int length = nucleotides.length();
		StringBuilder codon = new StringBuilder();
		for(int pos=frame; pos<length; pos++) {
			
			char nt = nucleotides.charAt(pos);
			
			// skip deletion
			if(nt == Constants.NULL.charAt(0)) {
				continue;
			}
			
			codon.append(nt);
			if(codon.length() == 3) {
				String codonStr = codon.toString();
				char aa = Codon.nuclToAmino(codonStr.toUpperCase());
				peptides.append(aa);
				codon.setLength(0);
			}
		}
		return peptides.toString();
	}
	
	public static String getReverseComplement (String nucleotide) {
		StringBuilder reverseComplementNTs = new StringBuilder(nucleotide);
		int length = nucleotide.length();
		for(int i=0; i<length; i++) {
			switch(reverseComplementNTs.charAt(i)) {
				case 'A': reverseComplementNTs.setCharAt(i, 'T'); break;
				case 'C': reverseComplementNTs.setCharAt(i, 'G'); break;
				case 'T': reverseComplementNTs.setCharAt(i, 'A'); break;
				case 'G': reverseComplementNTs.setCharAt(i, 'C'); break;
				default : break;
			}
		}
		
		return reverseComplementNTs.reverse().toString();
	}
	
}
