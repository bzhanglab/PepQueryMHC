package progistar.scan.run;

import java.io.IOException;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;

import progistar.scan.data.Constants;
import progistar.scan.data.Parameters;
import progistar.scan.function.CheckMemory;

public class Main {

	//-i test/benchmark/test.tsv -b test/benchmark/test.bam -m scan -s nucleotide -o test/benchmark/test.scan -@ 1
	//-i test/benchmark/test.tsv -b test/benchmark/C3N-02145.T.Aligned.sortedByCoord.out.bam -m scan -s nucleotide -o test/benchmark/test.scan -@ 4
	//-i test/benchmark/target_test.tsv -b test/benchmark/C3N-02145.T.Aligned.sortedByCoord.out.bam -m target -s peptide -o test/benchmark/target_test.scan -@ 4
	//-i test/benchmark/C3N_02145_T_nonreference.tsv -b test/benchmark/C3N-02145.T.Aligned.sortedByCoord.out.bam -m scan -s peptide -o test/benchmark/C3N_02145_T_nonreference.scan -@ 4
	
	//-i test/benchmark/C3N_02145_T_nonreference.tsv -b test/benchmark/C3N-02145.T.Aligned.sortedByCoord.out.bam -m target -s peptide -o test/benchmark/C3N_02145_T_nonreference.scan -@ 4
	
	
	public static void main(String[] args) throws IOException, InterruptedException {
		// performance metrics //
		long startTime = System.currentTimeMillis();
		/////////////////////////
		
		printDescription(args);
		parseModes(args);
		
		if(Parameters.mode.equalsIgnoreCase(Constants.MODE_SCAN) || Parameters.mode.equalsIgnoreCase(Constants.MODE_TARGET)) {
			MatchBAM.run(args);
		} else if(Parameters.mode.equalsIgnoreCase(Constants.MODE_ANNOTATE)) {
			Annotate.run(args);
		} else if(Parameters.mode.equalsIgnoreCase(Constants.MODE_FASTQ)) {
			MatchFASTQ.run(args);
		} else if(Parameters.mode.equalsIgnoreCase(Constants.MODE_EXTRACT)) {
			// TODO
			
		}
		
		// check peak memory
		Parameters.peakMemory = Math.max(Parameters.peakMemory, CheckMemory.checkUsedMemoryMB());
		
		long endTime = System.currentTimeMillis();
		System.out.println("Total Elapsed Time: "+(endTime-startTime)/1000+" sec");
		System.out.println("Estimated Peak Memory: "+Parameters.peakMemory +" MB");
	}
	
	
	/**
	 * 
	 * @param args
	 */
	public static void parseModes (String[] args) {
		CommandLine cmd = null;
		Options options = new Options();
		
		// select only interesting arguments
		String[] nArgs = new String[2];
		int nIdx = 0;
		for(int i=0; i<args.length; i++) {
			if(args[i].equalsIgnoreCase("-m") || args[i].equalsIgnoreCase("--mode")) {
				nArgs[nIdx] = args[i];
				nArgs[++nIdx] = args[++i];
			}
		}
		
		Option optionMode = Option.builder("m")
				.longOpt("mode").argName("scan|target|fastq|annotate|extract")
				.hasArg()
				.required(true)
				.desc("\"scan-mode\" counts all reads matching a given sequence by traversing all reads and annotates their genomic information. "
						+ "\n\"target-mode\" counts all reads matching a given sequence in a given genomic region."
						+ "\n\"fastq-mode\" counts all reads matching a given sequence by traversing all reads. "
						+ "\n\"annotate-mode\" annotates the best class code for a given genomic region (it requires GTF file)."
						+ "\n\"extract-mode\" extracts matched reads for each queried sequence and output them as a BAM file."
						+ "\n\"scan and target modes\" require .bai in advance.")
				.build();
		
		options.addOption(optionMode);
		
		CommandLineParser parser = new DefaultParser();
	    HelpFormatter helper = new HelpFormatter();
	    boolean isFail = false;
	    
		try {
		    cmd = parser.parse(options, nArgs);
		    
		    if(cmd.hasOption("m")) {
		    	Parameters.mode = cmd.getOptionValue("m");
		    	
		    	if( !(  Parameters.mode.equalsIgnoreCase(Constants.MODE_SCAN) || 
		    			Parameters.mode.equalsIgnoreCase(Constants.MODE_TARGET) ||
		    			Parameters.mode.equalsIgnoreCase(Constants.MODE_FASTQ) ||
		    			Parameters.mode.equalsIgnoreCase(Constants.MODE_ANNOTATE) ||
		    			Parameters.mode.equalsIgnoreCase(Constants.MODE_EXTRACT)) ) {
		    		isFail = true;
		    	}
		    }
		    
		} catch (Exception e) {
			isFail = true;
		}
		
		if(isFail) {
		    helper.printHelp("Usage:", options);
		    System.exit(0);
		}
	}
	
	public static void printDescription (String[] args) {
		System.out.println(Constants.NAME+" "+Constants.VERSION+" ("+Constants.RELEASE+")");
		System.out.println("running date: " + java.time.LocalDate.now());
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
