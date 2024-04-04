package progistar.scan.data;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;

import progistar.scan.function.Translator;

public class ParseRecord {

	private ParseRecord() {}

	/**
	 * id fr_sequence
	 * id ACGTGGAGT
	 * 
	 * 
	 * @param file
	 * @return
	 * @throws IOException
	 */
	public static ArrayList<Record> parse (File file) throws IOException {
		ArrayList<Record> records = new ArrayList<Record>();
		BufferedReader BR = new BufferedReader(new FileReader(file));
		String line = null;
		
		Record.header = BR.readLine(); // skip header
		
		String[] headerSplit = Record.header.split("\t");
		int obsSeqIdx = -1;
		int genomicLociCntIdx = -1;
		int genomicLociIdx = -1;
		for(int i=0; i<headerSplit.length; i++) {
			if(headerSplit[i].equalsIgnoreCase("ObservedNucleotide")) {
				obsSeqIdx = i;
			} else if(headerSplit[i].equalsIgnoreCase("GenomicLociCount")) {
				genomicLociCntIdx = i;
			} else if(headerSplit[i].equalsIgnoreCase("GenomicLoci")) {
				genomicLociIdx = i;
			} 
		}
		
		while((line = BR.readLine()) != null) {
			String[] fields = line.split("\t");
			
			String sequence = fields[obsSeqIdx];
			String genomicLociCnt = fields[genomicLociCntIdx];
			String genomicLoci = fields[genomicLociIdx];
			
			Record record = new Record();
			record.record = line;
			record.frSequence = sequence;
			record.genomicLociCount = Integer.parseInt(genomicLociCnt);
			record.rcSequence = Translator.getReverseComplement(sequence);
			
			String chr = "*";
			int start = -1;
			int end = -1;
			
			if(!genomicLoci.equalsIgnoreCase("-")) {
				String[] gLoc = genomicLoci.split("\\|");
				for(String gLocus : gLoc) {
					chr = gLocus.split("\\:")[0];
					if(start == -1) {
						start = Integer.parseInt(gLocus.split("\\:")[1].split("\\-")[0]);
					}
					end = Integer.parseInt(gLocus.split("\\:")[1].split("\\-")[1]);
				}
				
				if(chr.equalsIgnoreCase("chrx")) {
					chr = "chrX";
				} else if(chr.equalsIgnoreCase("chry")) {
					chr = "chrY";
				} else if(chr.equalsIgnoreCase("chrm")) {
					chr = "chrM";
				} else if(!chr.startsWith("chr")) {
					chr = chr.toUpperCase();
				}
			}
			
			record.chr = chr;
			record.start = start;
			record.end = end;
			
			records.add(record);
		}
		
		BR.close();
		return records;
	}
}
