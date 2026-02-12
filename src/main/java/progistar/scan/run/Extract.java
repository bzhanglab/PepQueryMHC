package progistar.scan.run;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import progistar.scan.data.BarcodeTable;
import progistar.scan.data.Codon;
import progistar.scan.data.Constants;
import progistar.scan.data.Parameters;
import progistar.scan.data.Phred;
import progistar.scan.data.SequenceRecord;
import progistar.scan.fileIO.ParseRecord;
import progistar.scan.function.BAMIndex;
import progistar.scan.function.CheckMemory;

public class Extract {

	public static void run(String[] args) throws IOException, InterruptedException {
		parseExtractModes(args);
		Codon.mapping();
		Phred.loadTable(); // load phred table
		// single cell barcode
		if(Parameters.isSingleCellMode) {
			BarcodeTable.load();
		}
		
		// check BAM index
		BAMIndex.index(Parameters.bamFile);
		
		ArrayList<SequenceRecord> records = ParseRecord.parse(Parameters.inputFile);
		
		//// Prepare tasks
		ArrayList<Task> tasks = new ArrayList<Task>();
		
		// auto strand detection
		if(Parameters.strandedness.equalsIgnoreCase(Constants.AUTO_STRANDED)) {
			tasks.addAll(Task.getStrandDetectionTask());
			ExecutorService executorService = Executors.newFixedThreadPool(Parameters.threadNum);
			List<Worker> callableExList = new ArrayList<>();
			for(int i=0; i<tasks.size(); i++) {
				Task task = tasks.get(i);
				callableExList.add(new Worker(task, tasks.size()));
			}
			
			// check peak memory
			Parameters.peakMemory = Math.max(Parameters.peakMemory, CheckMemory.checkUsedMemoryMB());
			
			// reset done count
			Worker.resetDoneCount();
			executorService.invokeAll(callableExList);
			executorService.shutdown();
			
			int R1F = 0;
			int R1R = 0;
			int R2F = 0;
			int R2R = 0;
			
			for(Task task :tasks) {
				R1F += task.R1F;
				R1R += task.R1R;
				R2F += task.R2F;
				R2R += task.R2R;
			}
			
			if(R1F*10 < R1R && R2F > R2R*10) {
				Parameters.strandedness = Constants.RF_STRANDED;
			} else if(R1F > R1R*10 && R2F*10 < R2R) {
				Parameters.strandedness = Constants.FR_STRANDED;
			}// for single end
			else if( (R1F + R2F) > 10 * (R1R + R2R) ) {
				Parameters.strandedness = Constants.F_STRANDED;
			} else if( 10 * (R1F + R2F) < (R1R + R2R) ) {
				Parameters.strandedness = Constants.R_STRANDED;
			} 
			else {
				Parameters.strandedness = Constants.NON_STRANDED;
			}
			
			System.out.println("Estimate strandedness");
			System.out.println("1F\t1R\t2F\t2R");
			System.out.println(R1F+"\t"+R1R+"\t"+R2F+"\t"+R2R);
			
			if(R1F+R1R+R2F+R2R == 0) {
				System.out.println("Fail to estimate stradedness!");
				System.out.println("It looks single-end RNA-seq experiement. Please specify strandedness.");
				System.exit(1);
			} else {
				System.out.println("Strandedness: "+Parameters.strandedness+"-stranded");
			}
			
			tasks.clear();
		}
		/////////////////////////////////////////////////////////////////
		
		// core algorithm
		// target mode
		Parameters.chunkSize = (records.size() / (10 * Parameters.threadNum) ) +1;
		tasks.addAll(Task.getExtractModeTasks(records, Parameters.chunkSize));
		//// sort tasks by descending order
		// Priority: Unmapped > Mapped
		Collections.sort(tasks);
		
		//// Enroll tasks on a thread pool
		ExecutorService executorService = Executors.newFixedThreadPool(Parameters.threadNum);
		List<Worker> callableExList = new ArrayList<>();
		for(int i=0; i<tasks.size(); i++) {
			Task task = tasks.get(i);
			callableExList.add(new Worker(task, tasks.size()));
		}
		
		// check peak memory
		Parameters.peakMemory = Math.max(Parameters.peakMemory, CheckMemory.checkUsedMemoryMB());
		
		// reset done count
		Worker.resetDoneCount();
		executorService.invokeAll(callableExList);
		executorService.shutdown();
		//// End of tasks
		
		System.out.println("Done all tasks!");
		// check peak memory
		Parameters.peakMemory = Math.max(Parameters.peakMemory, CheckMemory.checkUsedMemoryMB());

		// update peak memory from the tasks.
		for(Task task : tasks) {
			Parameters.peakMemory = Math.max(Parameters.peakMemory, task.peakMemory);
			
		}
		
		// write and delete bam files
		SamReader reader = SamReaderFactory.makeDefault().open(new File(Parameters.bamFile.getAbsolutePath()));
		SAMFileHeader samHeader = reader.getFileHeader();
		reader.close();
		
		// create BAM writer
		String outputExtractedBamPath = Parameters.outputBaseFilePath+"."+Parameters.mode+".bam";
        SAMFileWriter writer =
                new SAMFileWriterFactory()
                        .makeBAMWriter(samHeader, false, 
                        		new File(outputExtractedBamPath));
		
		for(Task task : tasks) {
			// for unmapped
			String[] filePaths = {task.getExtractedBamFileNameForMapped(),
					task.getExtractedBamFileNameForUnmapped()};
			for(String filePath : filePaths) {
				File tmpFile = new File(filePath);
				
				if(tmpFile.exists()) {
					reader = SamReaderFactory.makeDefault().open(tmpFile);
					for (SAMRecord record : reader) {
						writer.addAlignment(record);
					}
					reader.close();
					tmpFile.delete();
				}
			}
		}
		
		writer.close();
	}
	
	public static void parseExtractModes (String[] args) {
		
		// select only interesting arguments
		String[] nArgs = new String[args.length-2];
		int nIdx = 0;
		for(int i=0; i<args.length; i++) {
			if( args[i].equalsIgnoreCase("-i") || args[i].equalsIgnoreCase("--input") ||
				args[i].equalsIgnoreCase("-b") || args[i].equalsIgnoreCase("--bam") ||
				args[i].equalsIgnoreCase("-o") || args[i].equalsIgnoreCase("--output") ||
				args[i].equalsIgnoreCase("-@") || args[i].equalsIgnoreCase("--thread") ||
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
		
		Option optionThread = Option.builder("@")
				.longOpt("thread").argName("int")
				.hasArg()
				.required(false)
				.desc("the number of threads.")
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
		.addOption(optionThread)
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
		    
		    if(cmd.hasOption("@")) {
		    	Parameters.threadNum = Integer.parseInt(cmd.getOptionValue("@"));
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
			System.out.println("Threads: "+Parameters.threadNum);
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
