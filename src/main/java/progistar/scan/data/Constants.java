package progistar.scan.data;

public class Constants {
	public static final String NAME = "PepQueryMHC";
	public static final String VERSION = "v1.0.2";
	public static final String RELEASE = "April 29, 2025";
	
	public static final String MODE_TARGET = "target";
	public static final String MODE_SCAN = "scan";
	public static final String MODE_FASTQ = "fastq";
	public static final String MODE_ANNOTATE = "annotate";
	
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
	public static final int TYPE_FASTQ_MODE_TASK						= 5;
	
	// sequencing read file
	public static final byte SEQ_BAM			=	0b1;
	public static final byte SEQ_FASTQ_SINGLE	=	0b10;
	public static final byte SEQ_FASTQ_PAIRED	=	0b1100;
	
	//
	public static final int MAPPED_READS		=	1;
	public static final int UNMAPPED_READS		=	2;
	
	
	public static final String NULL				=	".";
	public static final byte NON_CODING_TRANSCRIPT	=	0;
	public static final byte CODING_TRANSCRIPT		=	1;
	
	public static final byte NCDS = 0;
	public static final byte CDS = 1;
	public static final byte UTR5 = 2;
	public static final byte UTR3 = 3;
	public static final byte INTRON = 4;
	public static final byte INTERGENIC = 5;
	
	public static final String MARK_PC = "PC";
	public static final String MARK_FS = "FS";
	public static final String MARK_NCRNA = "ncRNA";
	public static final String MARK_UTR5 = "5`-UTR";
	public static final String MARK_UTR3 = "3`-UTR";
	public static final String MARK_INTRON = "IR";
	public static final String MARK_ASRNA = "asRNA";
	public static final String MARK_INTERGENIC = "IGR";
	public static final String MARK_ES			=	"ES"; // for alternative splicing form
	public static final String MARK_UNKNOWN	= "Unknown";
	
	// Note that FRAME_X denots NO_FRAME.
	public static final byte FRAME_0		=	0;
	public static final byte FARME_1		=	1;
	public static final byte FRAME_2		=	2;
	public static final byte FRAME_X		=	3;
	
	// Penalty
	public static double PENALTY_ES					=	5;
	public static double PENALTY_5UTR				=	20;
	public static double PENALTY_3UTR				=	20;
	public static double PENALTY_FS					=	20;
	public static double PENALTY_ncRNA				=	20;
	public static double PENALTY_IR					=	30;
	public static double PENALTY_asRNA				=	50;
	public static double PENALTY_IGR				=	70;
	public static double PENALTY_UNMAP				=	100;
	public static double PENALTY_WARNING			=	Double.MAX_VALUE;
	
	public static final String DEFAULT_BARCODE_ID	=	"Undefined";
	public static final String NULL_BARCODE_ID		=	"Null";
	public static final String OTHER_BARCODE_ID		=	"Others";
	
	public static final byte SNP			=	0;
	public static final byte INS			=	1;
	public static final byte DEL			=	2;
	public static final byte CLP			=	3;
	
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
	 * F_STRANDED: forward - single-end
	 */
	public static final String F_STRANDED			=	"f";
	/**
	 * R_STRANDED: reverse - single-end
	 */
	public static final String R_STRANDED			=	"r";
	/**
	 * Determine strand-specific using XS and FLAGs
	 */
	public static final String AUTO_STRANDED			=	"auto";
	
	// Output fields for target and scan mode //
	public static final String MATCHED_LOCATION		= "Matched_location";
	public static final String MATCHED_MUTATIONS	= "Matched_mutations";
	public static final String MATCHED_STRAND		= "Matched_strand";
	public static final String MATCHED_PEPTIDE		= "Matched_peptide";
	public static final String MATCHED_NUCLEOTIDE	= "Matched_nucleotide";
	public static final String MATCHED_REFNUCLEOTIDE= "Matched_reference_nucleotide";
	public static final String MATCHED_READ_COUNT	= "Matched_read_count";
	public static final String MATCHED_RPHM			= "Matched_RPHM";
	public static final String MATCHED_NUM_LOCATION	= "Matched_num_locations";
	
	// Output fields for annotate //
	public static final String ANNOTATION_COUNT	= "Annotation_count";
	public static final String GENE_ID			= "Gene_id";
	public static final String GENE_NAME		= "Gene_name";
	public static final String GENE_STRAND		= "Gene_strand";
	public static final String GENE_TYPE		= "Gene_type";
	public static final String CLASS_CODE		= "Class_code";
	public static final String UNIQUE_CLASS_CODE= "Unique_class_code";
	public static final String WARNING_TAG		= "Warning_tag";
	
	// For annotation mode, no matched chr warning tag
	public static final String WARNING_TAG_NO_CHR	= "No_matched_reference_name";
	public static final String WARNING_TAG_NO_LOCATION	= "No_location_information";
	
	public static String getFullNameOfStrandedness (String strandedness) {
		String fullName = null;
		
		if(strandedness.equalsIgnoreCase(Constants.NON_STRANDED)) {
			fullName = "Non-stranded";
		} else if(strandedness.equalsIgnoreCase(Constants.F_STRANDED)) {
			fullName = "Forward-stranded (single-end)";
		} else if(strandedness.equalsIgnoreCase(Constants.R_STRANDED)) {
			fullName = "Reverse-stranded (single-end)";
		} else if(strandedness.equalsIgnoreCase(Constants.FR_STRANDED)) {
			fullName = "FR-stranded (paired-end)";
		} else if(strandedness.equalsIgnoreCase(Constants.RF_STRANDED)) {
			fullName = "RF-stranded (paired-end)";
		} else if(strandedness.equalsIgnoreCase(Constants.AUTO_STRANDED)) {
			fullName = "Auto-detection";
		}
		
		
		return fullName;
	}
	
}
