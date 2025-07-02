package progistar.scan.run;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
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

import htsjdk.samtools.fastq.FastqRecord;
import progistar.scan.data.BarcodeTable;
import progistar.scan.data.Codon;
import progistar.scan.data.Constants;
import progistar.scan.data.LibraryTable;
import progistar.scan.data.LocTable;
import progistar.scan.data.Parameters;
import progistar.scan.data.Phred;
import progistar.scan.data.SequenceRecord;
import progistar.scan.fileIO.ParseRecord;
import progistar.scan.fileIO.WriteOutput;
import progistar.scan.function.CheckMemory;

public class MatchFASTQ {

	public static void run(String[] args) throws IOException, InterruptedException {
		parseFASTQMode(args);
		Codon.mapping();
		Phred.loadTable(); // load phred table
		// single cell barcode
		if(Parameters.isSingleCellMode) {
			BarcodeTable.load();
		}
		
		// load library table
		if(Parameters.libFile != null) {
			LibraryTable.loadTable(Parameters.libFile);
		}
		
		ArrayList<SequenceRecord> records = ParseRecord.parse(Parameters.inputFile);
		ArrayList<Task> tasks = Task.getFASTQModeTasks(records);
		
		//// Enroll tasks on a thread pool
		int subAssistantThreadNum = 1;
		if(Parameters.sequencingFileType == Constants.SEQ_FASTQ_PAIRED) {
			subAssistantThreadNum = 2;
		}
		ExecutorService executorService = Executors.newFixedThreadPool(subAssistantThreadNum);
		List<Worker> callableExList = new ArrayList<>();
		for(int i=0; i<tasks.size(); i++) {
			Task task = tasks.get(i);
			callableExList.add(new Worker(task, tasks.size()));
		}
		// check peak memory
		Parameters.peakMemory = Math.max(Parameters.peakMemory, CheckMemory.checkUsedMemoryMB());
		executorService.invokeAll(callableExList);
		executorService.shutdown();
		
		// reset done count
		Worker.resetDoneCount();
		//// End of tasks
		System.out.println("Done all tasks!");
		// check peak memory
		Parameters.peakMemory = Math.max(Parameters.peakMemory, CheckMemory.checkUsedMemoryMB());

		// update peak memory from the tasks.
		for(Task task : tasks) {
			Parameters.peakMemory = Math.max(Parameters.peakMemory, task.peakMemory);
		}
		
		// calculate library size
		if(LibraryTable.isEmpty()) {
			for(Task task : tasks) {
				task.processedReads.forEach((barcode, count)->{
					Double libSize = LibraryTable.table.get(barcode);
					if(libSize == null) {
						libSize = .0;
					}
					libSize += count;
					LibraryTable.table.put(barcode, libSize);
				});
			}
		}
		
		// make location table
		LocTable locTable = new LocTable();
		
		// union information
		for(Task task : tasks) {
			task.locTable.table.forEach((sequence, lInfos) -> {
				lInfos.forEach((key, lInfo)->{
					locTable.putLocation(lInfo);
				});
			});
		}
		WriteOutput.writeMainOutput(records, Parameters.outputBaseFilePath, locTable);
		//WriteOutput.writeLocationLevelOutput(records, Parameters.outputBaseFilePath, locTable);
		WriteOutput.writePeptideLevelOutput(records, Parameters.outputBaseFilePath, locTable);
		
	}
	
	public static void count (ArrayList<FastqRecord> records) {
		for(FastqRecord record : records) {
			Task.allTrie.parseText(record.getReadString());
		}
	}

	public static void parseFASTQMode (String[] args) {
		// predefined
		Parameters.union = Constants.UNION_SUM;
		Parameters.count = Constants.COUNT_PRIMARY;
		
		// select only interesting arguments
		String[] nArgs = new String[args.length-2];
		int nIdx = 0;
		for(int i=0; i<args.length; i++) {
			if( args[i].equalsIgnoreCase("-i") || args[i].equalsIgnoreCase("--input") ||
				args[i].equalsIgnoreCase("-0") || args[i].equalsIgnoreCase("--fastq_single") ||
				args[i].equalsIgnoreCase("-1") || args[i].equalsIgnoreCase("--fastq_paired_1") ||
				args[i].equalsIgnoreCase("-2") || args[i].equalsIgnoreCase("--fastq_paired_2") ||
				args[i].equalsIgnoreCase("-o") || args[i].equalsIgnoreCase("--output") ||
				args[i].equalsIgnoreCase("-@") || args[i].equalsIgnoreCase("--thread") ||
				args[i].equalsIgnoreCase("-l") || args[i].equalsIgnoreCase("--lib_size") ||
				args[i].equalsIgnoreCase("-w") || args[i].equalsIgnoreCase("--white_list") ||
				args[i].equalsIgnoreCase("-p") || args[i].equalsIgnoreCase("--prob") ||
				args[i].equalsIgnoreCase("-u") || args[i].equalsIgnoreCase("--union") ||
				args[i].equalsIgnoreCase("-seq") || args[i].equalsIgnoreCase("--sequence_column_name") ||
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
		Option optionFastq0 = Option.builder("0")
				.longOpt("fastq_single").argName("fastq|fastq.gz")
				.hasArg()
				.required(false)
				.desc("First-end FASTQ file.")
				.build();
		
		Option optionFastq1 = Option.builder("1")
				.longOpt("fastq_paired_1").argName("fastq|fastq.gz")
				.hasArg()
				.required(false)
				.desc("First-end FASTQ file.")
				.build();
		
		Option optionFastq2 = Option.builder("2")
				.longOpt("fastq_paired_2").argName("fastq|fastq.gz")
				.hasArg()
				.required(false)
				.desc("Second-end FASTQ file.")
				.build();
		////////////////////////////////
		
		Option optionOutput = Option.builder("o")
				.longOpt("output").argName("file path")
				.hasArg()
				.required(true)
				.desc("output prefix path.")
				.build();
		

		Option optionStrandeness = Option.builder("s")
				.longOpt("strand").argName("non|fr|rf|f|r")
				.hasArg()
				.required(true)
				.desc("strand-specificity. non: non-stranded, fr: fr-second strand, rf: fr-first strand, f: forward strand for single-end, r: reverse strand for single-end.")
				.build();
		
		Option optionThread = Option.builder("@")
				.longOpt("thread").argName("int")
				.hasArg()
				.required(false)
				.desc("the number of threads.")
				.build();
		
		Option optionIL = Option.builder("e")
				.longOpt("equal").argName("")
				.required(false)
				.desc("consider that I is equivalent to L (only available in scan mode).")
				.build();
		
		Option optionLibSize = Option.builder("l")
				.longOpt("lib_size").argName("file path")
				.hasArg()
				.required(false)
				.desc("TSV file including library size information." +
						"\nIf this option is not specified, then it estimates the library size automatically. This estimation takes additional time for target mode.")
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
		
		Option optionSequenceColumnName = Option.builder("seq")
				.longOpt("sequence_column_name").argName("string")
				.hasArg()
				.required(false)
				.desc("specify sequence column name, case-insensitive (default is sequence).")
				.build();
		
		options.addOption(optionInput)
		.addOption(optionOutput)
		.addOption(optionStrandeness)
		.addOption(optionFastq0)
		.addOption(optionFastq1)
		.addOption(optionFastq2)
		.addOption(optionThread)
		.addOption(optionIL)
		.addOption(optionLibSize)
		.addOption(optionVerbose)
		.addOption(optionWhiteList)
		.addOption(optionSequenceColumnName)
		.addOption(optionROIThreshold);
		
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
		    if(cmd.hasOption("0")) {
		    	Parameters.fastq0File = new File(cmd.getOptionValue("0"));
		    	Parameters.sequencingFileType |= Constants.SEQ_FASTQ_SINGLE;
		    }
		    
		    if(cmd.hasOption("1")) {
		    	Parameters.fastq1File = new File(cmd.getOptionValue("1"));
		    	Parameters.sequencingFileType |= 0b100;
		    }
		    
		    if(cmd.hasOption("2")) {
		    	Parameters.fastq2File = new File(cmd.getOptionValue("2"));
		    	Parameters.sequencingFileType |= 0b1000;
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
		    
		    if(cmd.hasOption("u")) {
		    	// default is max.
		    	if(cmd.getOptionValue("u").equalsIgnoreCase("sum")) {
		    		Parameters.union = Constants.UNION_SUM;
		    	}
		    }
		    
		    if(cmd.hasOption("seq")) {
		    	Parameters.sequenceColumnName = cmd.getOptionValue("seq");
		    }
		    
		} catch (ParseException e) {
			System.out.println(e.getMessage());
			isFail = true;
		}
		
		// at least one sequencing read option should provide.
		if( !(	Parameters.sequencingFileType == Constants.SEQ_FASTQ_SINGLE || 
				Parameters.sequencingFileType == Constants.SEQ_FASTQ_PAIRED) ) {
			isFail = true;
		}
		
		if(isFail) {
		    helper.printHelp("Usage:", options);
		    System.exit(0);
		} else {
			System.out.println("Input file name: "+Parameters.inputFile.getAbsolutePath());
			
			if(Parameters.sequencingFileType == Constants.SEQ_FASTQ_SINGLE) {
				System.out.println("FASTQ file name: "+Parameters.fastq0File.getAbsolutePath());
			} else if(Parameters.sequencingFileType == Constants.SEQ_FASTQ_PAIRED) {
				System.out.println("FASTQ_1 file name: "+Parameters.fastq1File.getAbsolutePath());
				System.out.println("FASTQ_2 file name: "+Parameters.fastq2File.getAbsolutePath());
			}
			
			System.out.println("Output file name: "+Parameters.outputBaseFilePath);

			if(Parameters.whitelistFile != null) {
				System.out.println("White-list file name: "+Parameters.whitelistFile.getName() +" (single-cell mode)");
			}
			
			System.out.println("Sequence column name (case-insensitive): "+Parameters.sequenceColumnName);
			System.out.println("Strandedness: "+Constants.getFullNameOfStrandedness(Parameters.strandedness));
			System.out.println("Mode: "+Parameters.mode);
			System.out.println("Count: "+Parameters.count);
			System.out.println("Peptide level count: "+Parameters.union);
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
