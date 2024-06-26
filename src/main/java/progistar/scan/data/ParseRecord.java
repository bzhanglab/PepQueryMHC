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
		
		BAMSRecord.header = BR.readLine();
		BAMSRecord.fileName = Scan.bamFile.getName();
		
		String[] headerSplit = BAMSRecord.header.split("\t");
		int obsSeqIdx = -1;
		int genomicLociIdx = -1;
		int strandIdx = -1;
		
		for(int i=0; i<headerSplit.length; i++) {
			if(Scan.sequence.equalsIgnoreCase(Constants.SEQUENCE_NUCLEOTIDE) && 
					headerSplit[i].equalsIgnoreCase("sequence")) {
				obsSeqIdx = i;
			} else if(Scan.sequence.equalsIgnoreCase(Constants.SEQUENCE_PEPTIDE) && 
					headerSplit[i].equalsIgnoreCase("sequence")) {
				obsSeqIdx = i;
			} else if(headerSplit[i].equalsIgnoreCase("location")) {
				genomicLociIdx = i;
			} else if(headerSplit[i].equalsIgnoreCase("strand")) {
				strandIdx = i;
			} 
		}
		
		if(Scan.mode.equalsIgnoreCase(Constants.MODE_TARGET)) {

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
				
				String chr = "*"; // unmapped read mark in bam
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
		} else {
			while((line = BR.readLine()) != null) {
				String[] fields = line.split("\t");
				String sequence = fields[obsSeqIdx];
				
				BAMSRecord record = new BAMSRecord();
				record.sequence = sequence;
				record.strand = "*";
				record.location = "-";
				
				String key = record.getKey();
				
				BAMSRecord indexedRecord = indxedRecords.get(key);
				if(indexedRecord == null) {
					indexedRecord = record;
					indxedRecords.put(key, indexedRecord);
					records.add(indexedRecord);
				}
				indexedRecord.records.add(line);
			}
		}
		
		
		BR.close();
		return records;
	}
	
	/**
	 * For target mode
	 * 
	 * @param records
	 * @param file
	 * @throws IOException
	 */
	public static void writeRecords (ArrayList<BAMSRecord> records, File file) throws IOException {
		BufferedWriter BW = new BufferedWriter(new FileWriter(file));
		BufferedWriter BWPeptCount = new BufferedWriter(new FileWriter(file.getAbsolutePath()+".pept_count.tsv"));
		
		// write header
		BW.append(BAMSRecord.header+"\t"+BAMSRecord.fileName);
		BW.newLine();
		BWPeptCount.append("Sequence\tReadCount");
		BWPeptCount.newLine();
		
		Hashtable<String, Long> readCountsPeptLevel = new Hashtable<String, Long>();
		// write records
		for(int i=0; i<records.size(); i++) {
			BAMSRecord record = records.get(i);
			
			long readCnt = record.readCnt;
			for(int j=0; j<record.records.size(); j++) {
				// heavy version:
				BW.append(record.records.get(j)).append("\t"+readCnt);
				BW.newLine();
			}
			
			Long sumReads = readCountsPeptLevel.get(record.sequence);
			if(sumReads == null) {
				sumReads = 0L;
			}
			sumReads += readCnt;
			readCountsPeptLevel.put(record.sequence, sumReads);
		}
		
		readCountsPeptLevel.forEach((sequence, reads)->{
			try {
				BWPeptCount.append(sequence+"\t"+reads);
				BWPeptCount.newLine();
			}catch(IOException ioe) {
				
			}
 		});
		
		BW.close();
		BWPeptCount.close();
	}
	
	/**
	 * For scan mode
	 * 
	 * @param records
	 * @param file
	 * @param tasks
	 * @throws IOException
	 */
	public static void writeRecords (ArrayList<BAMSRecord> records, File file, ArrayList<Task> tasks) throws IOException {
		BufferedWriter BW = new BufferedWriter(new FileWriter(file));
		BufferedWriter BWGenomicTuple = new BufferedWriter(new FileWriter(file.getAbsolutePath()+".gloc.tsv"));
		BufferedWriter BWNotFound = new BufferedWriter(new FileWriter(file.getAbsolutePath()+".not_found.tsv"));
		BufferedWriter BWPeptCount = new BufferedWriter(new FileWriter(file.getAbsolutePath()+".pept_count.tsv"));
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
		BW.append(BAMSRecord.header+"\tLocation\tMutations\tStrand\tObsNucleotide\tObsPeptide\tRefNucleotide\tReadCount");
		BW.newLine();
		BWGenomicTuple.append("Sequence\tLocation\tStrand\tReadCount");
		BWGenomicTuple.newLine();
		
		BWNotFound.append(BAMSRecord.header+"\tLocation");
		BWNotFound.newLine();
		BWPeptCount.append("Sequence\tReadCount");
		BWPeptCount.newLine();
		
		// write records
		// unique observed sequence.
		Hashtable<String, Long> readCountsPeptLevel = new Hashtable<String, Long>();
		Hashtable<String, Long> readCountsTupleLevel = new Hashtable<String, Long>();
		Hashtable<String, Boolean> isUniqueCal = new Hashtable<String, Boolean>();
		
		for(int i=0; i<records.size(); i++) {
			BAMSRecord record = records.get(i);
			String sequence = record.sequence;
			
			// IL equal mode
			if(Scan.isILEqual) {
				sequence = sequence.replace("I", "L");
			}
			
			ArrayList<LocationInformation> locations = locTable.getLocations(sequence);
			
			// if we process with I=L option, 
			// then redundant peptides can be appeared and counted multiple times
			if(isUniqueCal.get(sequence) == null) {
				for(LocationInformation location : locations) {
					long readCount = location.readCount;
					// it must be calculated once!
					// peptide level count
					Long sumReads = readCountsPeptLevel.get(location.obsPeptide);
					if(sumReads == null) {
						sumReads = 0L;
					}
					sumReads += readCount;
					readCountsPeptLevel.put(location.obsPeptide, sumReads);
					
					// tuple level count
					String tupleKey = location.obsPeptide+"\t"+location.location+"\t"+location.strand;
					sumReads = readCountsTupleLevel.get(tupleKey);
					if(sumReads == null) {
						sumReads = 0L;
					}
					sumReads += readCount;
					readCountsTupleLevel.put(tupleKey, sumReads);
				}
				isUniqueCal.put(sequence, true);
			}
			
			// if there are duplicated records, then this size > 1
			// if there is no duplication, then this size = 1
			for(int j=0; j<record.records.size(); j++) {
				if(locations.size() == 0) {
					BWNotFound.append(record.records.get(j)).append("\tNot found");
					BWNotFound.newLine();
				} else {
					for(LocationInformation location : locations) {
						// full information (including genomic sequence)
						BW.append(record.records.get(j)).append("\t"+location.getRes());
						BW.newLine();
					}
				}
			}
			
		}
		BWNotFound.close();
		BW.close();
		
		readCountsPeptLevel.forEach((sequence, reads)->{
			try {
				BWPeptCount.append(sequence+"\t"+reads);
				BWPeptCount.newLine();
			}catch(IOException ioe) {
				
			}
 		});
		
		BWPeptCount.close();
		
		readCountsTupleLevel.forEach((tupleKey, reads)->{
			try {
				BWGenomicTuple.append(tupleKey+"\t"+reads);
				BWGenomicTuple.newLine();
			}catch(IOException ioe) {
				
			}
 		});
		BWGenomicTuple.close();
	}
}
