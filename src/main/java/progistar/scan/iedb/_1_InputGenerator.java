package progistar.scan.iedb;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;

public class _1_InputGenerator {

	public static void main(String[] args) throws IOException {
		File fileHLAI = new File("/Users/seunghyukchoi/Documents/1_Projects/2024_BamScan/IEDB_Epitope_Human_HLA_I.tsv");
		File fileHLAII = new File("/Users/seunghyukchoi/Documents/1_Projects/2024_BamScan/IEDB_Epitope_Human_HLA_II.tsv");
		
		BufferedReader BRHLAI = new BufferedReader(new FileReader(fileHLAI));
		BufferedReader BRHLAII = new BufferedReader(new FileReader(fileHLAII));
		BufferedWriter BW = new BufferedWriter(new FileWriter("IEDB.I+II.tsv"));
		String line = null;
		
		String header = BRHLAI.readLine();
		header = "Type\tPeptide\tLength\t"+header.replace("\"", "");
		BW.append(header);
		BW.newLine();
		while((line = BRHLAI.readLine()) != null) {
			String[] fields = line.split("\t");
			String peptide = fields[2].replace("\"", "");
			String sequence = peptide.split("\\s")[0];
			
			BW.append("HLA-I\t");
			BW.append(sequence).append("\t");
			BW.append(sequence.length()+"\t");
			BW.append(line.replace("\"", ""));
			BW.newLine();
		}
		
		BRHLAI.close();
		
		BRHLAII.readLine(); // skip header
		while((line = BRHLAII.readLine()) != null) {
			String[] fields = line.split("\t");
			String peptide = fields[2].replace("\"", "");
			String sequence = peptide.split("\\s")[0];
			
			BW.append("HLA-II\t");
			BW.append(sequence).append("\t");
			BW.append(sequence.length()+"\t");
			BW.append(line.replace("\"", ""));
			BW.newLine();
		}
		BRHLAII.close();
		
		BW.close();
	}
}
