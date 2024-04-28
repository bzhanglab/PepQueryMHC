package progistar.scan.data;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Hashtable;

import progistar.scan.function.Translator;
import progistar.scan.run.Scan;
import progistar.scan.run.Task;

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
		Hashtable<String, BAMSRecord> indxedRecords = new Hashtable<String, BAMSRecord>();
		String line = null;
		
		BAMSRecord.header = BR.readLine()+"\t"+Scan.bamFile.getName(); // skip header
		
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
			record.sequence = sequence;
			record.strand = strand;
			record.location = genomicLoci;
			
			if(Scan.sequence.equalsIgnoreCase(Constants.SEQUENCE_NUCLEOTIDE)) {
				if(strand.charAt(0) == '-') {
					// rc sequence of nucleotide
					record.sequence = Translator.getReverseComplement(record.sequence);
				}
			}
			
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
			
			String key = record.getKey();
			
			BAMSRecord indexedRecord = indxedRecords.get(key);
			if(indexedRecord == null) {
				indexedRecord = record;
				indxedRecords.put(key, indexedRecord);
				records.add(indexedRecord);
			}
			indexedRecord.records.add(line);
		}
		
		BR.close();
		return records;
	}
	
	public static void writeRecords (ArrayList<BAMSRecord> records, File file) throws IOException {
		BufferedWriter BW = new BufferedWriter(new FileWriter(file));
		
		// write header
		BW.append(BAMSRecord.header);
		BW.newLine();
		
		// write records
		for(int i=0; i<records.size(); i++) {
			BAMSRecord record = records.get(i);
			
			int readCnt = record.readCnt;
			for(int j=0; j<record.records.size(); j++) {
				// heavy version:
				BW.append(record.records.get(j)).append("\t"+readCnt);
				BW.newLine();
			}
		}
		
		BW.close();
	}
	
	public static void writeRecords (ArrayList<BAMSRecord> records, File file, ArrayList<Task> tasks) throws IOException {
		BufferedWriter BW = new BufferedWriter(new FileWriter(file));
		LocTable locTable = new LocTable();
		
		// union information
		for(Task task : tasks) {
			task.locTable.table.forEach((sequence, lInfos) -> {
				lInfos.forEach((key, lInfo)->{
					locTable.putLocation(lInfo);
				});
			});
		}
		
		// write header
		BW.append(BAMSRecord.header);
		BW.newLine();
		
		// write records
		for(int i=0; i<records.size(); i++) {
			BAMSRecord record = records.get(i);
			String sequence = record.sequence;
			int readCnt = record.readCnt;
			for(int j=0; j<record.records.size(); j++) {
				ArrayList<String> locations = locTable.getLocations(sequence);
				for(String location : locations) {
					BW.append(record.records.get(j)).append("\t"+location);
					BW.newLine();
				}
			}
		}
		
		BW.close();
	}
}
