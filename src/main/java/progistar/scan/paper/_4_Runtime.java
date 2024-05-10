package progistar.scan.paper;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Hashtable;

public class _4_Runtime {

	public static void main(String[] args) throws IOException {
		File[] files = new File("/Users/seunghyukchoi/Documents/1_Projects/2024_BamScan/LUAD_10samples").listFiles();
		Hashtable<String, Double> libSizes = new Hashtable<String, Double>();
		Hashtable<String, Hashtable<String, Hashtable<String, String>>> pHLAs = new Hashtable<String, Hashtable<String, Hashtable<String ,String>>>();
		
		for(File file : files) {
			if(file.getName().startsWith(".")) continue;
			if(!file.getName().endsWith(".libsize")) continue;
			
			BufferedReader BR = new BufferedReader(new FileReader(file));
			String bamName = file.getName().replace(".libsize", "");
			double value = Double.parseDouble(BR.readLine());
			
			libSizes.put(bamName, value);
			
			BR.close();
			
			BR = new BufferedReader(new FileReader(file.getAbsolutePath().replace(".bam.libsize", ".scan")));
			String line = null;
			
			BR.readLine(); // skip header
			while((line = BR.readLine()) != null) {
				String[] fields = line.split("\t");
				String peptide = fields[0];
				String location = fields[1];
				
				Hashtable<String, Hashtable<String, String>> peptides = pHLAs.get(bamName);
				if(peptides == null) {
					peptides = new Hashtable<String, Hashtable<String, String>>();
					pHLAs.put(bamName, peptides);
				}
				Hashtable<String, String> locations = peptides.get(peptide);
				if(locations == null) {
					locations = new Hashtable<String, String>();
					peptides.put(peptide, locations);
				}
				locations.put(location, "");
			}
			
			BR.close();
		}
		
		File logFile = new File("/Users/seunghyukchoi/Documents/1_Projects/2024_BamScan/LUAD_10samples/runtime.log");
		BufferedReader BR = new BufferedReader(new FileReader(logFile));
		String line = null;
		
		Hashtable<String, Double> runTimes = new Hashtable<String, Double>();
		String bamName = null;
		while((line = BR.readLine()) != null) {
			if(line.startsWith("run")) {
				bamName = line.split("\\s")[1];
			} else if(line.startsWith("Total")) {
				double value = Double.parseDouble(line.split("\\s")[3]);
				runTimes.put(bamName, value);
			}
		}
		
		BR.close();
		System.out.println("File\tLibSize\tTime\tpHLAs\tLocations");
		runTimes.forEach((bam, time)->{
			try {
				BufferedWriter BW = new BufferedWriter(new FileWriter(bam+".locs"));
				BW.append("Peptide\tLocations");
				BW.newLine();
				Double libSize = libSizes.get(bam);
				Hashtable<String, Hashtable<String, String>> peptides = pHLAs.get(bam);
				long[] nums = new long[2];
				nums[0] = peptides.size();
				peptides.forEach((peptide, locations)->{
					nums[1] += locations.size();
					try {
						BW.append(peptide+"\t"+locations.size());
						BW.newLine();
					}catch(IOException ioe) {
						
					}
				});
				
				System.out.println(bam+"\t"+libSize+"\t"+time+"\t"+nums[0]+"\t"+nums[1]);
				BW.close();
			}catch(IOException ioe) {
				
			}
			
		});
	}
}
