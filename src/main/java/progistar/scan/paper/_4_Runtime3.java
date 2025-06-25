package progistar.scan.paper;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Hashtable;

public class _4_Runtime3 {

	public static void main(String[] args) throws IOException {
		File log = new File("/Volumes/Papers/2024_BamScan/4.Performance_test/Stranded/2.BamQuery/10k_9mers_3/logs/BamQuery_Res_LUAD.log");
		
		BufferedReader BR = new BufferedReader(new FileReader(log));
		String line = null;
		
		Hashtable<String, Double> times = new Hashtable<String, Double>();
		
		double totalTime = 0;
		double sumOfCountTime = 0;
		System.out.println("Type\tTime");
		while((line = BR.readLine()) != null) {
			if(line.contains("Processed Bam File :")) {
				String[] fields = line.split("\\s");
				
				String sampleName = null;
				String min = null;
				for(int i=0; i<fields.length; i++) {
					if(fields[i].contains("star_salmon")) {
						sampleName = fields[i].replace(".", "").replace("_star_salmon", "");
					} else if(fields[i].contains(".")) {
						min = fields[i];
					}
				}
				
				times.put(sampleName, Double.parseDouble(min));
				System.out.println(sampleName+"\t"+min);
				sumOfCountTime += Double.parseDouble(min);
			} else if(line.contains("Total time run function BamQuery to end")){
				String[] fields = line.split("\\s");
				for(int i=0; i<fields.length; i++) {
					if(fields[i].contains(".")) {
						totalTime = Double.parseDouble(fields[i]);
					}
				}
			}
		}
		
		BR.close();
		System.out.println("Indexing\t"+(totalTime - sumOfCountTime));
		
	}
}
