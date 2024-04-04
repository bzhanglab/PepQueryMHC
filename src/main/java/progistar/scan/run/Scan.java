package progistar.scan.run;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.LinkedList;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;

import progistar.scan.data.BAMSRecord;
import progistar.scan.data.Constants;
import progistar.scan.data.ParseRecord;
import progistar.scan.data.Table;

public class Scan {

	public static String inputFilePath = null;
	public static String sheetFilePath   = null;
	public static String outputFilePath	   = null;
	public static String mode	=	Constants.MODE_TARGET;
	public static String sequence	=	Constants.SEQUENCE_PEPTIDE;

	public static int threadNum = 1;
	public static int chunkSize = 1000;
	
	
	public static void main(String[] args) throws IOException {
		long startTime = System.currentTimeMillis();
		parseOptions(args);
		
		ArrayList<BAMSRecord> records = ParseRecord.parse(new File(inputFilePath));
		
		LinkedList<Task> tasks = Task.divideTasks(records, chunkSize);
		Worker[] workers = new Worker[threadNum];
		
		int workerIdx = 1;
		while(!tasks.isEmpty()) {
			Task task = tasks.pollFirst();;
			boolean isAssigned = false;
			
			while(!isAssigned) {
				for(int j=0; j<workers.length; j++) {
					if(workers[j] == null || !workers[j].isAlive()) {
						workers[j] = new Worker(workerIdx++, task);
						workers[j].start();
						isAssigned = true;
						break;
					}
				}
				
				Thread.yield();
			}
		}
		
		boolean isProcessing = true;
		while(isProcessing) {
			isProcessing = false;
			for(int i=0; i<workers.length; i++) {
				if(workers[i] != null)
					isProcessing |= workers[i].isAlive();
			}
			Thread.yield();
		}
		
		
		System.out.println("Done all tasks!");
		Table table = new Table();

		for(Task task : tasks) {
			table.addStat(task);
		}
		table.write(new File(outputFilePath), tasks.get(0).records);
		
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
		
		Option optionOutput = Option.builder("o")
				.longOpt("output").argName("file path")
				.hasArg()
				.required(true)
				.desc("output matrix path")
				.build();
		
		Option optionMode = Option.builder("m")
				.longOpt("mode").argName("scan/target")
				.hasArg()
				.required(true)
				.desc("\"scan-mode\" counts all reads matching a given sequence by traversing all reads and annotates their genomic information. "
						+ "\n\"target-mode\" counts all reads matching a given sequence in a given genomic region."
						+ " \"target-mode\" requires .bai in advance.")
				.build();
		
		Option optionSequence = Option.builder("s")
				.longOpt("sequence").argName("peptide/nucleotide")
				.hasArg()
				.required(true)
				.desc("sequence type")
				.build();
		
		Option optionThread = Option.builder("@")
				.longOpt("thread").argName("int")
				.hasArg()
				.required(false)
				.desc("the number of threads")
				.build();
		
		
		options.addOption(optionInput)
		.addOption(optionOutput)
		.addOption(optionMode)
		.addOption(optionSequence)
		.addOption(optionThread);
		
		CommandLineParser parser = new DefaultParser();
	    HelpFormatter helper = new HelpFormatter();
	    boolean isFail = false;
	    
		try {
		    cmd = parser.parse(options, args);
		    
		    if(cmd.hasOption("i")) {
		    	inputFilePath = cmd.getOptionValue("i");
		    }
		    
		    if(cmd.hasOption("o")) {
		    	outputFilePath = cmd.getOptionValue("o");
		    }
		    
		    if(cmd.hasOption("m")) {
		    	mode = cmd.getOptionValue("m");
		    	
		    	if( !(mode.equalsIgnoreCase(Constants.MODE_FULL) || 
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
		    
		    if(cmd.hasOption("@")) {
		    	threadNum = Integer.parseInt(cmd.getOptionValue("@"));
		    }
		} catch (ParseException e) {
			System.out.println(e.getMessage());
			isFail = true;
		}
		
		if(isFail) {
		    helper.printHelp("Usage:", options);
		    System.exit(0);
		} else {
			System.out.println("Mode: "+mode);
			System.out.println("Threads: "+threadNum);
		}
		
		System.out.println();
	}
}
