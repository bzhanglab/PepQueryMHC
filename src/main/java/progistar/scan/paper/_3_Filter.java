package progistar.scan.paper;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Arrays;
import java.util.Hashtable;

public class _3_Filter {

	public static void main(String[] args) throws IOException {
		File file = new File("/Users/seunghyukchoi/Documents/1_Projects/2024_BamScan/IEDB.I+II.tsv");
		
		BufferedReader BR = new BufferedReader(new FileReader(file));
		BufferedWriter BW = new BufferedWriter(new FileWriter(file.getAbsolutePath().replace(".tsv",".lenfilter.repeatfilter.tsv")));
		String line = null;
		
		BW.append(BR.readLine());
		BW.newLine();
		
		Hashtable<String, String> removeRepeatSequence = new Hashtable<String, String>();
		int[] repeats = new int[26];
		while((line = BR.readLine()) != null) {
			String[] fields = line.split("\t");
			String peptide = fields[1];
			String hla = fields[0];
			int length = peptide.length();
			// repeat check
			Arrays.fill(repeats, 0);
			int diversity = 0;
			for(int i=0; i<peptide.length(); i++) {
				if(repeats[peptide.charAt(i)-'A'] == 0) {
					diversity ++;
					repeats[peptide.charAt(i)-'A'] = 1;
				}
			}
			if(diversity <= 2) {
				removeRepeatSequence.put(peptide, "");
				continue;
			}
			
			if(hla.equalsIgnoreCase("HLA-I") ) {
				if(length >= 8 && length <= 15) {
					BW.append(line);
					BW.newLine();
				}
			} else {
				if(length >= 8 && length <= 25) {
					BW.append(line);
					BW.newLine();
				}
			}
		}
		BW.close();
		BR.close();
		
		System.out.println(removeRepeatSequence.size()+" repeated sequences were removed");
	}
}
