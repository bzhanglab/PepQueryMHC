package progistar.scan.run;

import java.io.File;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;

import progistar.scan.data.Constants;
import progistar.scan.data.Parameters;

public class Extract {

	
	public static void parseExtractModes (String[] args) {
		
		// select only interesting arguments
		String[] nArgs = new String[args.length-2];
		int nIdx = 0;
		for(int i=0; i<args.length; i++) {
			if( args[i].equalsIgnoreCase("-i") || args[i].equalsIgnoreCase("--input") ||
				args[i].equalsIgnoreCase("-b") || args[i].equalsIgnoreCase("--bam") ||
				args[i].equalsIgnoreCase("-o") || args[i].equalsIgnoreCase("--output") ||
				args[i].equalsIgnoreCase("-c") || args[i].equalsIgnoreCase("--count") ||
				args[i].equalsIgnoreCase("-w") || args[i].equalsIgnoreCase("--white_list") ||
				args[i].equalsIgnoreCase("-p") || args[i].equalsIgnoreCase("--prob") ||
				args[i].equalsIgnoreCase("-seq") || args[i].equalsIgnoreCase("--sequence_column_name") ||
				args[i].equalsIgnoreCase("-loc") || args[i].equalsIgnoreCase("--location_column_name") ||
				args[i].equalsIgnoreCase("-s") || args[i].equalsIgnoreCase("--strand")) {
				nArgs[nIdx++] = args[i++];
				nArgs[nIdx++] = args[i];
			} 
			else if( args[i].equalsIgnoreCase("-v") || args[i].equalsIgnoreCase("--verbose") ||
					 args[i].equalsIgnoreCase("-e") || args[i].equalsIgnoreCase("--equal")) {
				nArgs[nIdx++] = args[i];
			}
		}
		
		
		CommandLine cmd = null;
		Options options = new Options();
		
		// Mandatory
		Option optionInput = Option.builder("i")
				.longOpt("input").argName("file path")
				.hasArg()
				.required(true)
				.desc("input path.")
				.build();
		
		//////////////// Sequencing reads /////////////////
		/**
		 * Sequencing reads must be provided!
		 * 
		 *
		 */
		Option optionBam = Option.builder("b")
				.longOpt("bam").argName("bam|sam")
				.hasArg()
				.required(true)
				.desc("bam or sam file.")
				.build();
		////////////////////////////////
		
		Option optionOutput = Option.builder("o")
				.longOpt("output").argName("file path")
				.hasArg()
				.required(true)
				.desc("output prefix path.")
				.build();
		
		Option optionPrimary = Option.builder("c")
				.longOpt("count").argName("primary|all")
				.hasArg()
				.required(false)
				.desc("count only primary or all reads (default is primary).")
				.build();
		
		Option optionIL = Option.builder("e")
				.longOpt("equal").argName("")
				.required(false)
				.desc("consider that I is equivalent to L (only available in scan mode).")
				.build();
		
		Option optionVerbose = Option.builder("v")
				.longOpt("verbose").argName("")
				.required(false)
				.desc("print every messages being processed.")
				.build();
		
		Option optionWhiteList = Option.builder("w")
				.longOpt("white_list").argName("file path")
				.hasArg()
				.required(false)
				.desc("cell barcode list (tsv).")
				.build();
		
		Option optionROIThreshold = Option.builder("p")
				.longOpt("prob").argName("float (0,1]")
				.hasArg()
				.required(false)
				.desc("ignore ROIs (region of interests) with greater than a given error probability (default is 0.05).")
				.build();
		
		Option optionStrandeness = Option.builder("s")
				.longOpt("strand").argName("non|fr|rf|f|r|auto")
				.hasArg()
				.required(false)
				.desc("strand-specificity. non: non-stranded, fr: fr-second strand, rf: fr-first strand, f: forward strand for single-end, r: reverse strand for single-end, "
						+ "auto: auto-detection. Auto-detection is only available if there is XS tag in a given BAM file (default is auto).")
				.build();
		
		Option optionSequenceColumnName = Option.builder("seq")
				.longOpt("sequence_column_name").argName("string")
				.hasArg()
				.required(false)
				.desc("specify sequence column name, case-insensitive (default is sequence).")
				.build();
		
		Option optionLocationColumnName = Option.builder("loc")
				.longOpt("location_column_name").argName("string")
				.hasArg()
				.required(false)
				.desc("specify location column name, case-insensitive (default is location).")
				.build();
		
		options.addOption(optionInput)
		.addOption(optionOutput)
		.addOption(optionStrandeness)
		.addOption(optionBam)
		.addOption(optionPrimary)
		.addOption(optionIL)
		.addOption(optionVerbose)
		.addOption(optionWhiteList)
		.addOption(optionROIThreshold)
		.addOption(optionSequenceColumnName)
		.addOption(optionLocationColumnName);
		
		CommandLineParser parser = new DefaultParser();
	    HelpFormatter helper = new HelpFormatter();
	    boolean isFail = false;
	    Parameters.sequencingFileType = 0b0;
	    
		try {
		    cmd = parser.parse(options, nArgs);
		    
		    if(cmd.hasOption("i")) {
		    	Parameters.inputFile = new File(cmd.getOptionValue("i"));
		    }
		    
		    ///////// input sequencing reads /////////
		    if(cmd.hasOption("b")) {
		    	Parameters.bamFile = new File(cmd.getOptionValue("b"));
		    	Parameters.sequencingFileType |= Constants.SEQ_BAM;
		    }
		    /////////////////////////////////////
		    
		    if(cmd.hasOption("o")) {
		    	Parameters.outputBaseFilePath = new File(cmd.getOptionValue("o")).getAbsolutePath();
		    }
		    
		    if(cmd.hasOption("s")) {
		    	Parameters.strandedness = cmd.getOptionValue("s");
		    	// there is no matched option
		    	if(!Parameters.strandedness.equalsIgnoreCase(Constants.AUTO_STRANDED) &&
		    		!Parameters.strandedness.equalsIgnoreCase(Constants.FR_STRANDED) &&
		    		!Parameters.strandedness.equalsIgnoreCase(Constants.RF_STRANDED) &&
		    		!Parameters.strandedness.equalsIgnoreCase(Constants.F_STRANDED) &&
		    		!Parameters.strandedness.equalsIgnoreCase(Constants.R_STRANDED) &&
		    		!Parameters.strandedness.equalsIgnoreCase(Constants.NON_STRANDED) ) {
		    		System.out.println("Wrong strandedness: "+Parameters.strandedness);
		    		isFail = true;
		    	}
		    }
		    
		    if(cmd.hasOption("e")) {
		    	Parameters.isILEqual = true;
		    }
		    
		    if(cmd.hasOption("c")) {
		    	Parameters.count = cmd.getOptionValue("c");
		    	
		    	if(Parameters.count.equalsIgnoreCase(Constants.COUNT_PRIMARY)) {
		    		Parameters.count = Constants.COUNT_PRIMARY;
		    	} else {
		    		Parameters.count = Constants.COUNT_ALL;
		    	}
		    }
		    
		    if(cmd.hasOption("v")) {
		    	Parameters.verbose = true;
		    }
		    
		    if(cmd.hasOption("l")) {
		    	Parameters.libFile = new File(cmd.getOptionValue("l"));
		    }
		    
		    if(cmd.hasOption("seq")) {
		    	Parameters.sequenceColumnName = cmd.getOptionValue("seq");
		    }
		    
		    if(cmd.hasOption("loc")) {
		    	Parameters.locationColumnName = cmd.getOptionValue("loc");
		    }
		    
		    if(cmd.hasOption("w")) {
		    	Parameters.whitelistFile = new File(cmd.getOptionValue("w"));
		    	Parameters.isSingleCellMode = true;
		    }
		    
		    if(cmd.hasOption("p")) {
		    	double roiCutoff = Double.parseDouble(cmd.getOptionValue("p"));
		    	if(Math.abs(roiCutoff) > 1 || roiCutoff == 0) {
		    		System.out.println("ROI cutoff is out of range (0,1]: "+roiCutoff);
		    		isFail = true;
		    	} else {
		    		Parameters.ROIErrorThreshold = roiCutoff;
		    	}
		    	
		    }
		    
		} catch (ParseException e) {
			System.out.println(e.getMessage());
			isFail = true;
		}
		
		if(isFail) {
		    helper.printHelp("Usage:", options);
		    System.exit(0);
		} else {
			System.out.println("Input file name: "+Parameters.inputFile.getAbsolutePath());
			System.out.println("BAM/SAM file name: "+Parameters.bamFile.getAbsolutePath());
			System.out.println("Output file name: "+Parameters.outputBaseFilePath);

			if(Parameters.whitelistFile != null) {
				System.out.println("White-list file name: "+Parameters.whitelistFile.getName() +" (single-cell mode)");
			}
			
			System.out.println("Sequence column name (case-insensitive): "+Parameters.sequenceColumnName);
			System.out.println("Location column name (case-insensitive): "+Parameters.locationColumnName);
			System.out.println("Strandedness: "+Constants.getFullNameOfStrandedness(Parameters.strandedness));
			System.out.println("Mode: "+Parameters.mode);
			System.out.println("Count: "+Parameters.count);
			System.out.println("ROI cutoff: "+Parameters.ROIErrorThreshold);
			if(Parameters.verbose) {
				System.out.println("Verbose messages");
			}
			if(Parameters.isILEqual) {
				if(Parameters.mode.equalsIgnoreCase(Constants.MODE_SCAN) && Parameters.sequence.equalsIgnoreCase(Constants.SEQUENCE_PEPTIDE)) {
					System.out.println("I and L are equivalent!");
				} else {
					System.out.println("This is target mode or nucleotide input. IL option is ignored.");
					Parameters.isILEqual = false;
				}
			}
		}
		System.out.println();
	}
}
