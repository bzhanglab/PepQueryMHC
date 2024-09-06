package progistar.scan.paper;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Hashtable;

public class _4_Runtime2 {

	public static void main(String[] args) throws IOException {
		File log = new File("/Volumes/Zhanglab/2024_BamScan/4.Performance_test/1.BamScan/prob_005/log");
		
		BufferedReader BR = new BufferedReader(new FileReader(log));
		String line = null;
		
		Hashtable<String, Double> times = new Hashtable<String, Double>();
		Hashtable<String, Double> mems = new Hashtable<String, Double>();
		
		String kmers = "";
		ArrayList<String> lists = new ArrayList<String>();
		while((line = BR.readLine()) != null) {
			if(line.startsWith("Records without")) {
				kmers = line.split("\\s")[line.split("\\s").length-1];
				if(times.get(kmers) == null) {
					lists.add(kmers);
				}
			} else if(line.startsWith("Total Elapsed Time")) {
				String time = line.split("\\s")[line.split("\\s").length-2];
				
				Double that = times.get(kmers);
				if(that == null) {
					that = .0;
				}
				
				that += Double.parseDouble(time);
				times.put(kmers, that);
			} else if(line.startsWith("Estimated Peak Memory")) {
				String mem = line.split("\\s")[line.split("\\s").length-2];
				
				Double that = mems.get(kmers);
				if(that == null) {
					that = .0;
				}
				
				that = Math.max(Double.parseDouble(mem), that);
				mems.put(kmers, that);
			}
		}
		
		BR.close();
		
		for(int i=0; i<lists.size(); i++) {
			kmers = lists.get(i);
			System.out.println(kmers+"\t"+times.get(kmers)+"\t"+mems.get(kmers));
			
		}
		
	}
}
