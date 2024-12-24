package progistar.scan.data;

import java.io.File;

import progistar.scan.function.CheckMemory;

public class Parameters {

	public static long peakMemory = CheckMemory.checkUsedMemoryMB();
	
	public static File inputFile = null;
	public static File bamFile = null;
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
}
