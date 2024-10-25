package progistar.scan.data;

public class Constants {
	public static final String NAME = "PepQueryHLA";
	public static final String VERSION = "v1.0.0b+";
	public static final String RELEASE = "September 24, 2024";
	
	public static final String MODE_TARGET = "target";
	public static final String MODE_SCAN = "scan";
	
	public static final String SEQUENCE_PEPTIDE = "peptide";
	/**
	 * @deprecated
	 */
	public static final String SEQUENCE_NUCLEOTIDE = "nucleotide";
	
	public static final String COUNT_ALL= "all";
	public static final String COUNT_PRIMARY = "primary";
	
	public static final String UNION_MAX = "max";
	public static final String UNION_SUM = "sum";
	
	public static final int	TYPE_TARGET_MODE_TASK						= 1;
	public static final int TYPE_TARGET_MODE_LIBRARY_ESTIMATION_TASK 	= 2;
	public static final int TYPE_SCAN_MODE_TASK							= 3;
	public static final int	TYPE_STRAND_DETECTION_TASK					= 4;
	
	//
	public static final int MAPPED_READS		=	1;
	public static final int UNMAPPED_READS		=	2;
	
	public static final byte EXON = 30;
	public static final byte CDS = 100;
	
	public static final String NULL				=	".";
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
	
	public static final byte SNP			=	0;
	public static final byte INS			=	1;
	public static final byte DEL			=	2;
	public static final byte CLP			=	3;
	
	public static final String DEFAULT_BARCODE_ID	=	"ReadCount";
	public static final String NULL_BARCODE_ID		=	"Null";
	public static final String OTHER_BARCODE_ID		=	"Others";
	
	// Strand-specific
	public static final String NON_STRANDED			=	"non";
	/**
	 * FR_STRANDED: fr-second strand, direct stranded
	 * the first read: forward
	 * the second read: reverse
	 */
	public static final String FR_STRANDED			=	"fr";
	/**
	 * RF_STRANDED: fr-first strand, reverse stranded
	 * the first read: reverse
	 * the second read: forward
	 */
	public static final String RF_STRANDED			=	"rf";
	/**
	 * F_STRANDED: forward - single-read
	 */
	public static final String F_STRANDED			=	"f";
	/**
	 * R_STRANDED: reverse - single-read
	 */
	public static final String R_STRANDED			=	"r";
	/**
	 * Determine strand-specific using XS and FLAGs
	 */
	public static final String AUTO_STRANDED			=	"auto";
	
	// Output fields //
	public static final String MATCHED_LOCATION		= "Matched_location";
	public static final String MATCHED_MUTATIONS	= "Matched_mutations";
	public static final String MATCHED_STRAND		= "Matched_strand";
	public static final String MATCHED_PEPTIDE		= "Matched_peptide";
	public static final String MATCHED_NUCLEOTIDE	= "Matched_nucleotide";
	public static final String MATCHED_REFNUCLEOTIDE= "Matched_reference_nucleotide";
	public static final String MATCHED_READ_COUNT	= "Matched_read_count";
	public static final String MATCHED_RPHM			= "Matched_RPHM";
	public static final String MATCHED_NUM_LOCATION	= "Matched_num_locations";
	
}
