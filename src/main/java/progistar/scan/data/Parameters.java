package progistar.scan.data;

import java.io.File;

import progistar.scan.function.CheckMemory;

public class Parameters {

	public static long peakMemory = CheckMemory.checkUsedMemoryMB();
	
	public static File inputFile = null;
	
	public static byte sequencingFileType = 0b0;
	// bam/sam file
	public static File bamFile = null;
	// for single-end fastq file
	public static File fastq0File = null;
	// for paired-ends fastq files
	public static File fastq1File = null;
	public static File fastq2File = null;
	
	public static File libFile = null;
	public static File whitelistFile = null;
	public static File gtfFile = null;
	
	public static String outputBaseFilePath	= null;
	public static String mode	=	Constants.MODE_TARGET;
	public static String sequence	=	Constants.SEQUENCE_PEPTIDE;
	public static String count	=	Constants.COUNT_PRIMARY;
	public static String union	=	Constants.UNION_SUM;
	public static String strandedness = Constants.AUTO_STRANDED;
	
	public static boolean isILEqual = false;
	public static boolean isSingleCellMode = false;
	public static boolean verbose = false;
	public static boolean stretch = false;
	public static boolean fullName = false;
	public static int threadNum = 4;
	public static int chunkSize = 100;
	
	
	public static int longestSequenceLen = -1;
	
	
	// read quality control ///////////////////
	/**
	 * single base cutoff
	 * @deprecated
	 */
	public static int singleBaseThreshold = 20;
	
	/**
	 * ROI base cutoff
	 */
	public static double ROIErrorThreshold = 0.05;
	///////////////////////////////////////////
	
	public static String unmmapedMarker = null;
	
	/////////////////// Annotate //////////////
	// this hinders assigning correct CDS.
	public static String[] drop_tag_list = {"cds_start_NF", "cds_end_NF"};
	
	// barcode separator in a read name
	public static String barcodeSeparatorInReadName = "_";
	
	
	public static String MARK_IF = Constants.MARK_SHORT_IF;
	public static String MARK_OOF = Constants.MARK_SHORT_OOF;
	public static String MARK_NCRNA = Constants.MARK_SHORT_NCRNA;
	public static String MARK_UTR5 = Constants.MARK_SHORT_UTR5;
	public static String MARK_UTR3 = Constants.MARK_SHORT_UTR3;
	public static String MARK_INTRON = Constants.MARK_SHORT_INTRON;
	public static String MARK_ASRNA = Constants.MARK_SHORT_ASRNA;
	public static String MARK_INTERGENIC = Constants.MARK_SHORT_INTERGENIC;
	public static String MARK_ES			=	Constants.MARK_SHORT_ES;
	public static String MARK_EE			=	Constants.MARK_SHORT_ASS;
	public static String MARK_UNKNOWN	= Constants.MARK_SHORT_UNKNOWN;
	
	public static void setFullName () {
		if(Parameters.fullName) {
			MARK_IF = Constants.MARK_FULL_IF;
			MARK_OOF = Constants.MARK_FULL_OOF;
			MARK_NCRNA = Constants.MARK_FULL_NCRNA;
			MARK_UTR5 = Constants.MARK_FULL_UTR5;
			MARK_UTR3 = Constants.MARK_FULL_UTR3;
			MARK_INTRON = Constants.MARK_FULL_INTRON;
			MARK_ASRNA = Constants.MARK_FULL_ASRNA;
			MARK_INTERGENIC = Constants.MARK_FULL_INTERGENIC;
			MARK_ES			=	Constants.MARK_FULL_ES;
			MARK_EE			=	Constants.MARK_FULL_EE;
			MARK_UNKNOWN	= Constants.MARK_FULL_UNKNOWN;
		}
 	}
	
	// input-columns
	public static String sequenceColumnName = "sequence";
	public static String locationColumnName = "location";
	public static String strandColumnName = "strand";
	
			
}
