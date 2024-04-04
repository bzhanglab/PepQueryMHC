package progistar.rna.data;

public class Constants {

	public static final String MODE_TARGET = "target";
	public static final String MODE_SCAN = "scan";
	
	public static final byte EXON = 30;
	public static final byte CDS = 100;
	
	public static final String ID_NULL				=	"-";
	public static final byte NON_CODING_TRANSCRIPT	=	0;
	public static final byte CODING_TRANSCRIPT		=	1;
	
	public static final byte UTR5 = 25;
	public static final byte UTR3 = 23;
	public static final byte NCDS = 20;
	public static final byte INTRON = 0;
	public static final byte INTERGENIC = -4;
	
	
	// Regional Character
	public static final char MARK_CDS 			=	'C';
	public static final char MARK_5UTR			=	'F';
	public static final char MARK_3UTR			=	'T';
	public static final char MARK_NCDS			=	'N';
	public static final char MARK_INTRON		=	'I';
	public static final char MARK_INTERGENIC	=	'-';
	
	
	// Note that FRAME_X denots NO_FRAME.
	public static final byte FRAME_0		=	0;
	public static final byte FARME_1		=	1;
	public static final byte FRAME_2		=	2;
	public static final byte FRAME_X		=	3;
	
	// Alternative Splicing Character
	public static final char MARK_AS			=	'A'; // for alternative splicing form
	public static final char MARK_CA			=	'C'; // for canonical form
}
