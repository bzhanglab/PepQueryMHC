package progistar.rna.run;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;

public class MakeSheet {
	
	public static String inputFilePath = null;
	public static String outputFilePath = null;
	public static String removeStr = null;
	public static String threadNum = "4";

	public static void main(String[] args) throws IOException, InterruptedException {
		parseOptions(args);
		
		BufferedWriter BW = new BufferedWriter(new FileWriter(outputFilePath));
		File[] files = new File(inputFilePath).listFiles();
		
		BW.append("sample_name\tfile_path\ttotal_read_count");
		BW.newLine();
		for(File file :files) {
			if(file.getName().startsWith(".")) continue;
			if(file.getName().endsWith(".bam")) {
				System.out.println(file.getName()+" is processed...");
				String sampleName = file.getName().replace(removeStr, "");
				String[] command = { "samtools", "view", "-F", "256", "-c", file.getAbsolutePath(), "-@", threadNum };
				
				ProcessBuilder processBuilder = new ProcessBuilder(command);

		        // Redirect error stream to output stream
		        processBuilder.redirectErrorStream(true);

		        // Start the process
		        Process process = processBuilder.start();

		        // Read output
		        BufferedReader reader = new BufferedReader(new InputStreamReader(process.getInputStream()));
		        String line;
		        StringBuilder output = new StringBuilder();
		        while ((line = reader.readLine()) != null) {
		            output.append(line);
		        }

		        // Wait for the process to finish and get the exit code
		        int exitCode = process.waitFor();
				long totalReadCount= Long.parseLong(output.toString());
				BW.append(sampleName+"\t"+file.getAbsolutePath()+"\t"+totalReadCount);
				BW.newLine();
			}
		}
		
		BW.close();
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
		
		Option optionSuffix = Option.builder("r")
				.longOpt("remove").argName("string")
				.hasArg()
				.required(true)
				.desc("sample name is written as a file name by removing a given string")
				.build();
		
		Option optionThread = Option.builder("@")
				.longOpt("thread").argName("number")
				.hasArg()
				.required(false)
				.desc("the number of threads to process. Default value is 4.")
				.build();
		
		
		options.addOption(optionInput)
		.addOption(optionOutput)
		.addOption(optionThread)
		.addOption(optionSuffix);
		
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
		    
		    if(cmd.hasOption("r")) {
		    	removeStr = cmd.getOptionValue("r");
		    }
		    
		    if(cmd.hasOption("@")) {
		    	threadNum = cmd.getOptionValue("@");
		    }
		} catch (ParseException e) {
			System.out.println(e.getMessage());
			isFail = true;
		}
		
		if(isFail) {
		    helper.printHelp("Usage:", options);
		    System.exit(0);
		} else {
		}
		
	}
}
