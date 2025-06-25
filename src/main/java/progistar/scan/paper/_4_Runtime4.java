package progistar.scan.paper;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Hashtable;

public class _4_Runtime4 {

	public static void main(String[] args) throws IOException {
		File log = new File("/Volumes/Papers/2024_BamScan/4.Performance_test/Stranded/1.PepQueryHLA/run.log");
		
		BufferedReader BR = new BufferedReader(new FileReader(log));
		String line = null;
		
		String replicate = ".10k.3";
		System.out.println("Type\tTime");
		
		String sampleName = null;
		while((line = BR.readLine()) != null) {
			if(line.contains("[OUT]")) {
				String[] fields = line.split("\\s");
				sampleName = fields[1];
			} else if(line.contains("Total Elapsed Time:")) {
				String[] fields = line.split("\\s");
				double time = Double.parseDouble(fields[3]);
				
				if(sampleName.contains(replicate)) {
					sampleName = sampleName.replace(replicate, "");
					System.out.println(sampleName+"\t"+(time/60));
				}
			}
		}
		
		BR.close();
		
	}
}
