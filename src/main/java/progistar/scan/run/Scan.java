package progistar.scan.run;

import java.io.File;
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

import progistar.scan.data.Codon;
import progistar.scan.data.Constants;
import progistar.scan.data.ParseRecord;
import progistar.scan.data.SequenceRecord;

public class Scan {

	//-i test/benchmark/test.tsv -b test/benchmark/test.bam -m scan -s nucleotide -o test/benchmark/test.scan -@ 1
	//-i test/benchmark/test.tsv -b test/benchmark/C3N-02145.T.Aligned.sortedByCoord.out.bam -m scan -s nucleotide -o test/benchmark/test.scan -@ 4
	//-i test/benchmark/target_test.tsv -b test/benchmark/C3N-02145.T.Aligned.sortedByCoord.out.bam -m target -s peptide -o test/benchmark/target_test.scan -@ 4
	//-i test/benchmark/C3N_02145_T_nonreference.tsv -b test/benchmark/C3N-02145.T.Aligned.sortedByCoord.out.bam -m scan -s peptide -o test/benchmark/C3N_02145_T_nonreference.scan -@ 4
	
	//-i test/benchmark/C3N_02145_T_nonreference.tsv -b test/benchmark/C3N-02145.T.Aligned.sortedByCoord.out.bam -m target -s peptide -o test/benchmark/C3N_02145_T_nonreference.scan -@ 4
	
	
	
	public static File inputFile = null;
	public static File bamFile = null;
	public static File outputFile	   = null;
	public static String mode	=	Constants.MODE_TARGET;
	public static String sequence	=	Constants.SEQUENCE_PEPTIDE;
	public static String count	=	Constants.COUNT_PRIMARY;
	public static double libSize = 0;
	public static boolean isILEqual = false;

	public static boolean verbose = false;
	public static boolean isRandom = false;
	public static int threadNum = 4;
	public static int chunkSize = 100;
	public static String unmmapedMarker = null;
	
	public static void main(String[] args) throws IOException, InterruptedException {
		
		long startTime = System.currentTimeMillis();
		printDescription(args);
		parseOptions(args);
		Codon.mapping();
		
		ArrayList<SequenceRecord> records = ParseRecord.parse(inputFile);
		
		//// Prepare tasks
		ArrayList<Task> tasks = new ArrayList<Task>();
		if(mode.equalsIgnoreCase(Constants.MODE_TARGET)) {
			// estimate library size
			if(libSize == 0) {
				tasks.addAll(Task.getLibSizeTask(records));
			}
			// target mode
			chunkSize = (records.size() / (10 * threadNum) ) +1;
			tasks.addAll(Task.getTargetModeTasks(records, chunkSize));
		} else if(mode.equalsIgnoreCase(Constants.MODE_SCAN)) {
			tasks.addAll(Task.getScanModeTasks(records));
		}
		//// sort tasks by descending order
		// Priority: Library > Unmapped > Mapped
		Collections.sort(tasks);
		
		//// Enroll tasks on a thread pool
		ExecutorService executorService = Executors.newFixedThreadPool(threadNum);
		List<Worker> callableExList = new ArrayList<>();
		for(int i=0; i<tasks.size(); i++) {
			Task task = tasks.get(i);
			callableExList.add(new Worker(task));
		}
		
		executorService.invokeAll(callableExList);
		executorService.shutdown();
		//// End of tasks
		
		System.out.println("Done all tasks!");
		
		// calculate library size
		if(libSize == 0) {
			for(Task task : tasks) {
				libSize += task.processedReads;
			}
		}
		
		if(mode.equalsIgnoreCase(Constants.MODE_TARGET)) {
			ParseRecord.writeRecords(records, outputFile);
		} else if(mode.equalsIgnoreCase(Constants.MODE_SCAN)) {
			ParseRecord.writeRecords(records, outputFile, tasks);
		}
		
		long endTime = System.currentTimeMillis();
		System.out.println("Total Elapsed Time: "+(endTime-startTime)/1000+" sec");
	}
	
	
	
	/**
	 * Parse and apply arguments
	 * 
	 * @param args
	 */
	public static void parseOptions (String[] args) {
		
		CommandLine cmd = null;
		Options options = new Options();
		
		// Mandatory
		Option optionInput = Option.builder("i")
				.longOpt("input").argName("file path")
				.hasArg()
				.required(true)
				.desc("input path")
				.build();
		
		Option optionBam = Option.builder("b")
				.longOpt("bam").argName("bam/sam")
				.hasArg()
				.required(true)
				.desc("bam or sam file")
				.build();
		
		Option optionOutput = Option.builder("o")
				.longOpt("output").argName("file path")
				.hasArg()
				.required(true)
				.desc("output matrix path")
				.build();
		
		Option optionMode = Option.builder("m")
				.longOpt("mode").argName("scan/target")
				.hasArg()
				.required(false)
				.desc("\"scan-mode\" counts all reads matching a given sequence by traversing all reads and annotates their genomic information. "
						+ "\n\"target-mode\" counts all reads matching a given sequence in a given genomic region."
						+ " \"target-mode\" requires .bai in advance.")
				.build();
		
		Option optionSequence = Option.builder("s")
				.longOpt("sequence").argName("peptide/nucleotide")
				.hasArg()
				.required(false)
				.desc("sequence type")
				.build();
		
		Option optionThread = Option.builder("@")
				.longOpt("thread").argName("int")
				.hasArg()
				.required(false)
				.desc("the number of threads")
				.build();
		
		Option optionPrimary = Option.builder("c")
				.longOpt("count").argName("primary/all")
				.hasArg()
				.required(false)
				.desc("count only primary or all reads")
				.build();
		
		Option optionIL = Option.builder("e")
				.longOpt("equal").argName("")
				.required(false)
				.desc("consider that I is equivalent to L (only available in scan mode).")
				.build();
		
		Option optionLibSize = Option.builder("l")
				.longOpt("lib_size").argName("int")
				.hasArg()
				.required(false)
				.desc("library size to calculate RPHM value." +
						"\nif this option is not used, then it estimates the library size automatically. This estimation takes additional time for target mode.")
				.build();
		
		Option optionRandomDist = Option.builder("r")
				.longOpt("random").argName("")
				.required(false)
				.desc("generate and match reversed sequence to find a random distribution (only available in scan mode).")
				.build();
		
		Option optionVerbose = Option.builder("v")
				.longOpt("verbose").argName("")
				.required(false)
				.desc("print every messages being processed.")
				.build();
		
		options.addOption(optionInput)
		.addOption(optionOutput)
		.addOption(optionMode)
		.addOption(optionSequence)
		.addOption(optionBam)
		.addOption(optionThread)
		.addOption(optionPrimary)
		.addOption(optionIL)
		.addOption(optionLibSize)
		.addOption(optionRandomDist)
		.addOption(optionVerbose);
		
		CommandLineParser parser = new DefaultParser();
	    HelpFormatter helper = new HelpFormatter();
	    boolean isFail = false;
	    
		try {
		    cmd = parser.parse(options, args);
		    
		    if(cmd.hasOption("i")) {
		    	inputFile = new File(cmd.getOptionValue("i"));
		    }
		    
		    if(cmd.hasOption("b")) {
		    	bamFile = new File(cmd.getOptionValue("b"));
		    }
		    
		    if(cmd.hasOption("o")) {
		    	outputFile = new File(cmd.getOptionValue("o"));
		    }
		    
		    if(cmd.hasOption("m")) {
		    	mode = cmd.getOptionValue("m");
		    	
		    	if( !(mode.equalsIgnoreCase(Constants.MODE_SCAN) || 
		    	   mode.equalsIgnoreCase(Constants.MODE_TARGET)) ) {
		    		isFail = true;
		    	}
		    }
		    
		    if(cmd.hasOption("s")) {
		    	sequence = cmd.getOptionValue("s");
		    	
		    	if( !(sequence.equalsIgnoreCase(Constants.SEQUENCE_NUCLEOTIDE) || 
		    			sequence.equalsIgnoreCase(Constants.SEQUENCE_PEPTIDE)) ) {
		    		isFail = true;
		    	}
		    }
		    
		    if(cmd.hasOption("e")) {
		    	isILEqual = true;
		    }
		    
		    if(cmd.hasOption("r")) {
		    	isRandom = true;
		    }
		    
		    if(cmd.hasOption("@")) {
		    	threadNum = Integer.parseInt(cmd.getOptionValue("@"));
		    }
		    
		    if(cmd.hasOption("c")) {
		    	count = cmd.getOptionValue("c");
		    	
		    	if(count.equalsIgnoreCase(Constants.COUNT_PRIMARY)) {
		    		count = Constants.COUNT_PRIMARY;
		    	} else {
		    		count = Constants.COUNT_ALL;
		    	}
		    }
		    
		    if(cmd.hasOption("v")) {
		    	verbose = true;
		    }
		    
		    if(cmd.hasOption("l")) {
		    	libSize = Double.parseDouble(cmd.getOptionValue("l"));
		    }
		    
		} catch (ParseException e) {
			System.out.println(e.getMessage());
			isFail = true;
		}
		
		if(isFail) {
		    helper.printHelp("Usage:", options);
		    System.exit(0);
		} else {
			System.out.println("Input file name: "+inputFile.getName());
			System.out.println("BAM/SAM file name: "+bamFile.getName());
			System.out.println("Output file name: "+outputFile.getName());
			System.out.println("Type: "+sequence);
			System.out.println("Mode: "+mode);
			System.out.println("Count: "+count);
			System.out.println("Threads: "+threadNum);
			if(verbose) {
				System.out.println("Verbose messages");
			}
			if(isILEqual) {
				if(mode.equalsIgnoreCase(Constants.MODE_SCAN) && sequence.equalsIgnoreCase(Constants.SEQUENCE_PEPTIDE)) {
					System.out.println("I and L are equivalent!");
				} else {
					System.out.println("This is target mode or nucleotide input. IL option is ignored.");
					isILEqual = false;
				}
			}
			
			if(isRandom) {
				if(mode.equalsIgnoreCase(Constants.MODE_SCAN) && sequence.equalsIgnoreCase(Constants.SEQUENCE_PEPTIDE)) {
					System.out.println("Generate random sequences.");
				} else {
					System.out.println("This is target mode or nucleotide input. Generation of random distribution is ignored.");
					isRandom = false;
				}
			}
		}
		System.out.println();
	}
	
	public static void printDescription (String[] args) {
		System.out.println(Constants.NAME+" "+Constants.VERSION+" (running date: " + java.time.LocalDate.now()+")");
		StringBuilder optionStr = new StringBuilder();
		optionStr.append("command line: ");
		for(int i=0; i<args.length; i++) {
			if(i != 0) {
				optionStr.append(" ");
			}
			optionStr.append(args[i]);
		}
		System.out.println(optionStr.toString());
		System.out.println();
	}
}
