package progistar.scan.data;

public class Constants {
	public static final String NAME = "PepQueryMHC";
	public static final String VERSION = "v1.0.2a";
	public static final String RELEASE = "May 6, 2025";
	
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
	
	public static final String MARK_FULL_IF = "In-frame";
	public static final String MARK_FULL_OOF = "Out-of-frame";
	public static final String MARK_FULL_NCRNA = "Non-coding RNA";
	public static final String MARK_FULL_UTR5 = "5`-untranslated region";
	public static final String MARK_FULL_UTR3 = "3`-untranslated region";
	public static final String MARK_FULL_INTRON = "intron retention";
	public static final String MARK_FULL_ASRNA = "Antisense RNA";
	public static final String MARK_FULL_INTERGENIC = "Intergenic region";
	public static final String MARK_FULL_ES			=	"Exon-skipping"; // for exon-skipping
	public static final String MARK_FULL_EE			=	"Alternative-splicing-site"; // for ASS
	public static final String MARK_FULL_UNKNOWN	= "Unknown";
	
	public static final String MARK_SHORT_IF = "IF";
	public static final String MARK_SHORT_OOF = "OOF";
	public static final String MARK_SHORT_NCRNA = "ncRNA";
	public static final String MARK_SHORT_UTR5 = "5`-UTR";
	public static final String MARK_SHORT_UTR3 = "3`-UTR";
	public static final String MARK_SHORT_INTRON = "IR";
	public static final String MARK_SHORT_ASRNA = "asRNA";
	public static final String MARK_SHORT_INTERGENIC = "IGR";
	public static final String MARK_SHORT_ES			=	"ES"; // for exon-skipping
	public static final String MARK_SHORT_EE			=	"ASS"; // for exon-exclusion
	public static final String MARK_SHORT_UNKNOWN	= "Unknown";
	
	// Note that FRAME_X denots NO_FRAME.
	public static final byte FRAME_0		=	0;
	public static final byte FARME_1		=	1;
	public static final byte FRAME_2		=	2;
	public static final byte FRAME_X		=	3;
	
	// Penalty
	//// Use between transcripts
	public static int PENALTY_ES				=	15;
	public static int PENALTY_EE				=	15;
	public static int PENALTY_5UTR				=	20;
	public static int PENALTY_3UTR				=	20;
	public static int PENALTY_OOF				=	20;
	public static int PENALTY_NCRNA				=	30;
	public static int PENALTY_IR				=	60;
	public static int PENALTY_ASRNA				=	120;
	public static int PENALTY_IGR				=	240;
	public static int PENALTY_UNMAP				=	480;
	public static int PENALTY_WARNING			=	1000;
	
	// Annotation priority: Select lower values as a representative.
	/// Use to present annotation for each transcript
	//// For regional variations
	public static int ANNOTATION_PRIORITY_IF		=	7;
	public static int ANNOTATION_PRIORITY_OOF		=	6;
	public static int ANNOTATION_PRIORITY_5UTR		=	5;
	public static int ANNOTATION_PRIORITY_3UTR		=	5;
	public static int ANNOTATION_PRIORITY_NCRNA		=	5;
	public static int ANNOTATION_PRIORITY_IR		=	4;
	public static int ANNOTATION_PRIORITY_ASRNA		=	3;
	public static int ANNOTATION_PRIORITY_IGR		=	2;
	public static int ANNOTATION_PRIORITY_UNMAP		=	1;
	//// For structural variations
	public static int ANNOTATION_PRIORITY_ES	=	1;
	public static int ANNOTATION_PRIORITY_EE	=	1;
	
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
