package progistar.scan.paper;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.Hashtable;

public class _4_Runtime {

	public static void main(String[] args) throws IOException {
		File[] files = new File("/Volumes/Zhanglab/2024_BamScan/4.Performance_test/1.BamScan").listFiles();
		Hashtable<String, Double> libSizes = new Hashtable<String, Double>();
		
		for(File file : files) {
			if(file.getName().startsWith(".")) continue;
			if(!file.getName().contains("0.05")) continue;
			if(!file.getName().endsWith(".libsize")) continue;
			
			BufferedReader BR = new BufferedReader(new FileReader(file));
			String bamName = file.getName().replace(".libsize", "");
			double value = Double.parseDouble(BR.readLine());
			
			libSizes.put(bamName, value);
			
			BR.close();
		}
		
		File logFile = new File("/Volumes/Zhanglab/2024_BamScan/4.Performance_test/1.BamScan/log");
		BufferedReader BR = new BufferedReader(new FileReader(logFile));
		String line = null;
		
		Hashtable<String, Double> runTimes = new Hashtable<String, Double>();
		Hashtable<String, Double> peakMemory = new Hashtable<String, Double>();
		String bamName = null;
		while((line = BR.readLine()) != null) {
			if(line.startsWith("command")) {
				bamName = line.split("\\s")[line.split("\\s").length-1];
			} else if(line.startsWith("Total")) {
				double value = Double.parseDouble(line.split("\\s")[3]);
				runTimes.put(bamName, value);
			} else if(line.startsWith("Estimated Peak Memory: ")) {
				double value = Double.parseDouble(line.split("\\s")[3]);
				peakMemory.put(bamName, value);
			}
		}
		
		BR.close();
		System.out.println("File\tLibSize\tTime\tpHLAs\tLocations");
		runTimes.forEach((bam, time)->{
			Double libSize = libSizes.get(bam);
			Double peakMem = peakMemory.get(bam);
			if(libSize != null) {
				System.out.println(bam+"\t"+libSize+"\t"+time+"\t"+peakMem);
			}
			
		});
	}
}
