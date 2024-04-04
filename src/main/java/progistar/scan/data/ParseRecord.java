package progistar.scan.data;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;

import progistar.scan.run.Scan;

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
	public static ArrayList<BAMSRecord> parse (File file) throws IOException {
		ArrayList<BAMSRecord> records = new ArrayList<BAMSRecord>();
		BufferedReader BR = new BufferedReader(new FileReader(file));
		String line = null;
		
		BAMSRecord.header = BR.readLine(); // skip header
		
		String[] headerSplit = BAMSRecord.header.split("\t");
		int obsSeqIdx = -1;
		int genomicLociIdx = -1;
		int strandIdx = -1;
		
		for(int i=0; i<headerSplit.length; i++) {
			if(Scan.sequence.equalsIgnoreCase(Constants.SEQUENCE_NUCLEOTIDE) && 
					headerSplit[i].equalsIgnoreCase("nucleotide")) {
				obsSeqIdx = i;
			} else if(Scan.sequence.equalsIgnoreCase(Constants.SEQUENCE_PEPTIDE) && 
					headerSplit[i].equalsIgnoreCase("peptide")) {
				obsSeqIdx = i;
			} else if(headerSplit[i].equalsIgnoreCase("location")) {
				genomicLociIdx = i;
			} else if(headerSplit[i].equalsIgnoreCase("strand")) {
				strandIdx = i;
			} 
		}
		
		while((line = BR.readLine()) != null) {
			String[] fields = line.split("\t");
			
			String sequence = fields[obsSeqIdx];
			String genomicLoci = fields[genomicLociIdx];
			String strand = fields[strandIdx];
			
			BAMSRecord record = new BAMSRecord();
			record.record = line;
			record.sequence = sequence;
			record.strand = strand;
			
			String chr = "*";
			int start = -1;
			int end = -1;
			
			if(!genomicLoci.equalsIgnoreCase("*")) {
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
