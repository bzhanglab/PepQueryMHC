package progistar.scan.fileIO;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Hashtable;
import java.util.Iterator;

import progistar.scan.data.BarcodeTable;
import progistar.scan.data.Constants;
import progistar.scan.data.LibraryTable;
import progistar.scan.data.LocTable;
import progistar.scan.data.LocationInformation;
import progistar.scan.data.Parameters;
import progistar.scan.data.SequenceRecord;
import progistar.scan.function.Utils;
import progistar.scan.function.Validation;
import progistar.scan.function.WriteStatistics;

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
	public static ArrayList<SequenceRecord> parse (File file) throws IOException {
		ArrayList<SequenceRecord> records = new ArrayList<SequenceRecord>();
		BufferedReader BR = new BufferedReader(new FileReader(file));
		Hashtable<String, SequenceRecord> indexedRecords = new Hashtable<String, SequenceRecord>();
		String line = null;
		
		SequenceRecord.header = BR.readLine();
		
		if(Parameters.mode.equalsIgnoreCase(Constants.MODE_TARGET) || Parameters.mode.equalsIgnoreCase(Constants.MODE_SCAN)) {
			SequenceRecord.fileName = Parameters.bamFile.getName();
		}
		
		String[] headerSplit = SequenceRecord.header.split("\t");
		int inputSeqIdx = -1;
		int genomicLociIdx = -1;
		int strandIdx = -1;
		
		for(int i=0; i<headerSplit.length; i++) {
			if(Parameters.sequence.equalsIgnoreCase(Constants.SEQUENCE_NUCLEOTIDE) && 
					headerSplit[i].equalsIgnoreCase("sequence")) {
				inputSeqIdx = i;
			} else if(Parameters.sequence.equalsIgnoreCase(Constants.SEQUENCE_PEPTIDE) && 
					headerSplit[i].equalsIgnoreCase("sequence")) {
				inputSeqIdx = i;
			} else if(headerSplit[i].equalsIgnoreCase("location")) {
				genomicLociIdx = i;
			} else if(headerSplit[i].equalsIgnoreCase("strand")) {
				strandIdx = i;
			} 
		}
		
		
		// fail to find fields
		if(Parameters.mode.equalsIgnoreCase(Constants.MODE_TARGET)) {
			if(inputSeqIdx == -1 || genomicLociIdx == -1 || strandIdx == -1) {
				System.out.println("Fail to find column names: sequence, location, strand");
				System.out.println("Please specify exact column names in your input file.");
				System.exit(1);
			}
			
		} else if(Parameters.mode.equalsIgnoreCase(Constants.MODE_ANNOTATE)) {
			if(genomicLociIdx == -1 || strandIdx == -1) {
				System.out.println("Fail to find column names: location, strand");
				System.out.println("Please specify exact column names in your input file.");
				System.exit(1);
			}
		} else if(Parameters.mode.equalsIgnoreCase(Constants.MODE_SCAN)) {
			if(inputSeqIdx == -1) {
				System.out.println("Fail to find column names: sequence");
				System.out.println("Please specify exact column names in your input file.");
				System.exit(1);
			}
		}
			
		
		if(Parameters.mode.equalsIgnoreCase(Constants.MODE_TARGET) || Parameters.mode.equalsIgnoreCase(Constants.MODE_ANNOTATE)) {

			while((line = BR.readLine()) != null) {
				String[] fields = line.split("\t");
				String genomicLoci = fields[genomicLociIdx];
				String strand = fields[strandIdx];
				
				
				SequenceRecord record = new SequenceRecord();
				
				if(Parameters.mode.equalsIgnoreCase(Constants.MODE_TARGET)) {
					// if the sequence is invalid
					String sequence = fields[inputSeqIdx];
					if(Validation.checkValidSequence(sequence)) {
						record.sequence = sequence;
					} else {
						continue;
					}
				} else if(Parameters.mode.equalsIgnoreCase(Constants.MODE_ANNOTATE)) {
					// pass
				}
				
				record.strand = strand;
				record.location = genomicLoci;
				
				/** @deprecated
				if(Scan.sequence.equalsIgnoreCase(Constants.SEQUENCE_NUCLEOTIDE)) {
					if(strand.charAt(0) == '-') {
						// rc sequence of nucleotide
						record.sequence = Translator.getReverseComplement(record.sequence);
					}
				}
				**/
				
				String chr = Constants.NULL;
				int start = -1;
				int end = -1;
				
				if(!genomicLoci.equalsIgnoreCase(Constants.NULL)) {
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
				
				SequenceRecord indexedRecord = indexedRecords.get(key);
				if(indexedRecord == null) {
					indexedRecord = record;
					indexedRecords.put(key, indexedRecord);
					records.add(indexedRecord);
				}
				indexedRecord.records.add(line);
			}
		} else if (Parameters.mode.equalsIgnoreCase(Constants.MODE_SCAN) || Parameters.mode.equalsIgnoreCase(Constants.MODE_FASTQ)){
			while((line = BR.readLine()) != null) {
				String[] fields = line.split("\t");
				String sequence = fields[inputSeqIdx];
				
				// if the sequence is invalid
				if(!Validation.checkValidSequence(sequence)) {
					continue;
				}
				
				SequenceRecord record = new SequenceRecord();
				record.sequence = Parameters.isILEqual ? sequence.replace("I", "L") : sequence;
				record.strand = Constants.NULL;
				record.location = Constants.NULL;
				
				String key = record.getKey();
				
				SequenceRecord indexedRecord = indexedRecords.get(key);
				if(indexedRecord == null) {
					indexedRecord = record;
					indexedRecords.put(key, indexedRecord);
					records.add(indexedRecord);
				}
				indexedRecord.records.add(line);
			}
		}
		
		if(Parameters.mode.equalsIgnoreCase(Constants.MODE_TARGET) || 
			Parameters.mode.equalsIgnoreCase(Constants.MODE_SCAN) ||
			Parameters.mode.equalsIgnoreCase(Constants.MODE_FASTQ)){
			// set longestLengthOfInputSequences
			for(SequenceRecord record : records) {
				Parameters.longestSequenceLen = Math.max(record.sequence.length(), Parameters.longestSequenceLen);
			}
			System.out.println("Longest length of input sequences: "+Parameters.longestSequenceLen);
		} else if(Parameters.mode.equalsIgnoreCase(Constants.MODE_ANNOTATE)) {
			
		}
		System.out.println("Records without duplication: "+records.size());
		
		
		BR.close();
		return records;
	}
	
}
