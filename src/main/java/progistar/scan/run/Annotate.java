package progistar.scan.run;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Hashtable;
import java.util.LinkedList;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;

import progistar.scan.data.Annotation;
import progistar.scan.data.Constants;
import progistar.scan.data.Gene;
import progistar.scan.data.GeneArray;
import progistar.scan.data.Parameters;
import progistar.scan.data.SequenceRecord;
import progistar.scan.fileIO.ParseGTF;
import progistar.scan.fileIO.ParseRecord;
import progistar.scan.fileIO.WriteOutput;
import progistar.scan.function.CheckMemory;
import progistar.scan.function.IndexConvertor;

public class Annotate {

	public static void run(String[] args) throws IOException {
		parseAnnotate(args);
		
		ArrayList<SequenceRecord> records = ParseRecord.parse(Parameters.inputFile);
		GeneArray[] geneArrays = ParseGTF.parseGTF(Parameters.gtfFile);
		
		Parameters.peakMemory = Math.max(Parameters.peakMemory, CheckMemory.checkUsedMemoryMB());
		
		Hashtable<String, LinkedList<Annotation>> allAnnotations = new Hashtable<String, LinkedList<Annotation>>();
		for(SequenceRecord sRecord : records) {
			LinkedList<Annotation> annotations = new LinkedList<Annotation>();
			int chrIdx = IndexConvertor.chrToIndex(sRecord.chr);
			
			// If the location is ".", then it is considered as "unknown"
			// or
			// no chr index found
			if(sRecord.location.equalsIgnoreCase(Constants.NULL) || chrIdx == -1) {
				// unmapped reads, unidentified region (.)
				Annotation annotation = new Annotation();
				annotation.classCode = Constants.MARK_UNKNOWN;
				
				// set warning tag
				if(sRecord.location.equalsIgnoreCase(Constants.NULL)) {
					annotation.warningTag = Constants.WARNING_TAG_NO_LOCATION;
				} else if(chrIdx == -1) {
					annotation.warningTag = Constants.WARNING_TAG_NO_CHR;
				}
				
				annotations.add(annotation);
			} 
			// location is given
			else {
				ArrayList<Gene> genes = geneArrays[chrIdx].findOverlap(sRecord.start, sRecord.end);
				if(genes.size() == 0) {
					// intergenic
					Annotation annotation = new Annotation();
					annotation.classCode = Constants.MARK_INTERGENIC;
					annotations.add(annotation);
				} else {
					for(Gene gene : genes) {
						Annotation[] annotates = gene.annotate(sRecord);
						for(Annotation annotation : annotates) {
							annotations.add(annotation);
						}
					}
				}
			}
			
			// calculate penalty
			annotations.forEach((annotation)->{
				annotation.calPenalty();
			});
			
			Collections.sort(annotations);
			
			// select annotations with the lowest penalty
			double minPenalty = annotations.peekFirst().penalty;
			while(true) {
				Annotation annotation = annotations.peekLast();
				if(annotation.penalty == minPenalty) {
					break;
				} else {
					// remove
					annotations.pollLast();
				}
			}
			
			// gene level summation
			annotations = Annotation.removeRedundancy(annotations);
			
			allAnnotations.put(sRecord.getKey(), annotations);
			Parameters.peakMemory = Math.max(Parameters.peakMemory, CheckMemory.checkUsedMemoryMB());
		}
		
		WriteOutput.writeAnnotateOutput(records, allAnnotations);
	}
	
	

	/**
	 * Parse and apply arguments
	 * 
	 * @param args
	 */
	public static void parseAnnotate (String[] args) {
		
		// select only interesting arguments
		String[] nArgs = new String[args.length-2];
		int nIdx = 0;
		for(int i=0; i<args.length; i++) {
			if( args[i].equalsIgnoreCase("-i") || args[i].equalsIgnoreCase("--input") ||
				args[i].equalsIgnoreCase("-g") || args[i].equalsIgnoreCase("--gtf") ||
				args[i].equalsIgnoreCase("-o") || args[i].equalsIgnoreCase("--output")) {
				nArgs[nIdx++] = args[i++];
				nArgs[nIdx++] = args[i];
			} 
			else if( args[i].equalsIgnoreCase("-v") || args[i].equalsIgnoreCase("--verbose") ||
					 args[i].equalsIgnoreCase("-s") || args[i].equalsIgnoreCase("--stretch")) {
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
		
		Option optionOutput = Option.builder("o")
				.longOpt("output").argName("file path")
				.hasArg()
				.required(true)
				.desc("output prefix path.")
				.build();
		
		Option optionGTF = Option.builder("g")
				.longOpt("gtf").argName("file path")
				.hasArg()
				.required(true)
				.desc("specify a GTF file to annotate a given genomic region.")
				.build();
		
		Option optionStretch = Option.builder("s")
				.longOpt("stretch").argName("")
				.required(false)
				.desc("enable to output single line per annotation.")
				.build();
		
		Option optionVerbose = Option.builder("v")
				.longOpt("verbose").argName("")
				.required(false)
				.desc("print every messages being processed.")
				.build();
		
		options.addOption(optionInput)
		.addOption(optionOutput)
		.addOption(optionVerbose)
		.addOption(optionStretch)
		.addOption(optionGTF);
		
		CommandLineParser parser = new DefaultParser();
	    HelpFormatter helper = new HelpFormatter();
	    boolean isFail = false;
	    
		try {
		    cmd = parser.parse(options, nArgs);
		    
		    if(cmd.hasOption("i")) {
		    	Parameters.inputFile = new File(cmd.getOptionValue("i"));
		    }
		    if(cmd.hasOption("o")) {
		    	Parameters.outputBaseFilePath = new File(cmd.getOptionValue("o")).getAbsolutePath();
		    }
		    
		    if(cmd.hasOption("g")) {
		    	Parameters.gtfFile = new File(cmd.getOptionValue("g"));
		    }
		    
		    if(cmd.hasOption("s")) {
		    	Parameters.stretch = true;
		    }
		    
		    if(cmd.hasOption("v")) {
		    	Parameters.verbose = true;
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
			System.out.println("Output file name: "+Parameters.outputBaseFilePath);
			System.out.println("GTF file name: "+Parameters.gtfFile.getAbsolutePath());

			System.out.println("Mode: "+Parameters.mode);
			if(Parameters.stretch) {
				System.out.println("Stretch output");
			}
			if(Parameters.verbose) {
				System.out.println("Verbose messages");
			}
		}
		System.out.println();
	}
	
}
