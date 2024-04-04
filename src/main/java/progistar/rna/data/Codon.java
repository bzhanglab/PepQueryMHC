package progistar.rna.data;

import java.util.Arrays;

public class Codon {
	private static final int nucleoIndexes = 8;
	private static String AminoToNuclArray[][];
	private static char NuclToAminoArray[][][];
	private static boolean setOkay = false;
	private static char aminoAcids[] = {'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y'};
	
	private static String nucleotides[][] = {
			/*A*/ {"GCT", "GCC", "GCA", "GCG"},
					{},
			/*C*/ {"TGT", "TGC"},
			/*D*/ {"GAT", "GAC"},
			/*E*/ {"GAA", "GAG"},
			/*F*/ {"TTT", "TTC"},
			/*G*/ {"GGT", "GGC", "GGA", "GGG"},
			/*H*/ {"CAT", "CAC"},
			/*I*/ {"ATT", "ATC", "ATA"},
					{},
			/*K*/ {"AAA", "AAG"},
			/*L*/ {"TTA", "TTG", "CTT", "CTC", "CTA", "CTG"},
			/*M*/ {"ATG"},
			/*N*/ {"AAT", "AAC"},
					{},
			/*P*/ {"CCT", "CCC", "CCA", "CCG"},
			/*Q*/ {"CAA", "CAG"},
			/*R*/ {"CGT", "CGC", "CGA", "CGG", "AGA", "AGG"},
			/*S*/ {"TCT", "TCC", "TCA", "TCG", "AGT", "AGC"},
			/*T*/ {"ACT", "ACC", "ACA", "ACG"},
					{},
			/*V*/ {"GTT", "GTA", "GTC", "GTG"},
			/*W*/ {"TGG"},
					{},
			/*Y*/ {"TAT", "TAC"},
					{}
	};
	
	
	public static void mapping() {
		AminoToNuclArray = new String[26][];
		NuclToAminoArray = new char[nucleoIndexes][nucleoIndexes][nucleoIndexes];
		
		for(int ntPos = 0; ntPos<nucleoIndexes; ntPos++){
			for(int ntPos_ = 0; ntPos_<nucleoIndexes; ntPos_++){
				Arrays.fill(NuclToAminoArray[ntPos][ntPos_], 'X');
			}
		}
		
		for(char AA : aminoAcids){
			AminoToNuclArray[AA - 'A'] = new String[nucleotides[AA -'A'].length];
			for(int ntPos = 0; ntPos<nucleotides[AA -'A'].length; ntPos++){
				AminoToNuclArray[AA - 'A'][ntPos] = nucleotides[AA -'A'][ntPos];
				
				NuclToAminoArray
				[nucleotides[AA -'A'][ntPos].charAt(0) & 7]
				[nucleotides[AA -'A'][ntPos].charAt(1) & 7]
				[nucleotides[AA -'A'][ntPos].charAt(2) & 7]
						= AA;
			}
		}
		
		setOkay = true;
	}
	
	
	public static Character nuclToAmino (String nucleotides){
		if(!setOkay) mapping();
		if(nucleotides.length() != 3) return 'X';
		return NuclToAminoArray[nucleotides.charAt(0) & 7][nucleotides.charAt(1) & 7][nucleotides.charAt(2) & 7];
	}
}