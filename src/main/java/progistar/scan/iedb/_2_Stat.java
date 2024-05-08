package progistar.scan.iedb;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.Hashtable;

public class _2_Stat {

	public static void main(String[] args) throws IOException {
		File file = new File("/Users/seunghyukchoi/Documents/1_Projects/2024_BamScan/IEDB.I+II.tsv");
		
		BufferedReader BR = new BufferedReader(new FileReader(file));
		String line = null;
		
		Hashtable<String, String> HLAIs = new Hashtable<String, String>();
		Hashtable<String, String> HLAIIs = new Hashtable<String, String>();
		int maxLength = 0;
		while((line = BR.readLine()) != null) {
			String[] fields = line.split("\t");
			String hla = fields[0];
			
			if(fields[1].length() == 2) {
				System.out.println(line);
			}
			
			if(hla.equalsIgnoreCase("HLA-I")) {
				HLAIs.put(fields[1], "");
			} else {
				HLAIIs.put(fields[1], "");
			}
			
			maxLength = Math.max(maxLength, fields[1].length());
		}
		
		BR.close();
		
		int[] HLAICounts = new int[maxLength+1];
		int[] HLAIICounts = new int[maxLength+1];
		
		HLAIs.forEach((sequence,nil)->{
			HLAICounts[sequence.length()]++;
		});
		HLAIIs.forEach((sequence,nil)->{
			HLAIICounts[sequence.length()]++;
		});
		
		for(int i=0; i<maxLength+1; i++) {
			if(i==0) {
				System.out.println("Length\tHLA-I\tHLA-II");
			} else {
				int hlaI = HLAICounts[i];
				int hlaII = HLAIICounts[i];
				System.out.println(i+"\t"+hlaI+"\t"+hlaII);
			}
		}
	}
}








