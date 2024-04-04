package progistar.rna.run;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;

import progistar.rna.data.Constants;
import progistar.rna.data.ParseRecord;
import progistar.rna.data.Sheet;
import progistar.rna.data.Table;

public class Scan {

	public static String inputFilePath = null;
	public static String sheetFilePath   = null;
	public static String outputFilePath	   = null;
	public static String mode	=	Constants.MODE_TARGET;
	public static int threadNum = 8;
	
	
	public static void main(String[] args) throws IOException {
		long startTime = System.currentTimeMillis();
		parseOptions(args);
		
		Sheet sheet = new Sheet(new File(sheetFilePath));
		
		ArrayList<Task> tasks = new ArrayList<Task>();
		Worker[] workers = new Worker[threadNum];
		
		int workerIdx = 1;
		for(int i=0; i<sheet.fileInfos.size(); i++) {
			
			System.out.println("Ready "+sheet.fileInfos.get(i).file.getName());
			Task task = new Task(ParseRecord.parse(new File(inputFilePath)), sheet.fileInfos.get(i));
			tasks.add(task);
			
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
		
		Option optionSheet = Option.builder("s")
				.longOpt("sheet").argName("tsv")
				.hasArg()
				.required(true)
				.desc("Sample information file.")
				.build();
		
		Option optionThread = Option.builder("@")
				.longOpt("thread").argName("int")
				.hasArg()
				.required(false)
				.desc("the number of threads")
				.build();
		
		
		options.addOption(optionInput)
		.addOption(optionOutput)
		.addOption(optionSheet)
		.addOption(optionMode)
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
		    
		    if(cmd.hasOption("s")) {
		    	sheetFilePath = cmd.getOptionValue("s");
		    }
		    
		    if(cmd.hasOption("m")) {
		    	mode = cmd.getOptionValue("m");
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
