package progistar.scan.data;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Hashtable;
import java.util.Iterator;

import progistar.scan.function.Random;
import progistar.scan.function.Translator;
import progistar.scan.function.Utils;
import progistar.scan.function.Validation;
import progistar.scan.function.WriteStatistics;
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
	public static ArrayList<SequenceRecord> parse (File file) throws IOException {
		ArrayList<SequenceRecord> records = new ArrayList<SequenceRecord>();
		BufferedReader BR = new BufferedReader(new FileReader(file));
		Hashtable<String, SequenceRecord> indexedRecords = new Hashtable<String, SequenceRecord>();
		String line = null;
		
		SequenceRecord.header = BR.readLine();
		SequenceRecord.fileName = Scan.bamFile.getName();
		
		String[] headerSplit = SequenceRecord.header.split("\t");
		int inputSeqIdx = -1;
		int genomicLociIdx = -1;
		int strandIdx = -1;
		
		for(int i=0; i<headerSplit.length; i++) {
			if(Scan.sequence.equalsIgnoreCase(Constants.SEQUENCE_NUCLEOTIDE) && 
					headerSplit[i].equalsIgnoreCase("sequence")) {
				inputSeqIdx = i;
			} else if(Scan.sequence.equalsIgnoreCase(Constants.SEQUENCE_PEPTIDE) && 
					headerSplit[i].equalsIgnoreCase("sequence")) {
				inputSeqIdx = i;
			} else if(headerSplit[i].equalsIgnoreCase("location")) {
				genomicLociIdx = i;
			} else if(headerSplit[i].equalsIgnoreCase("strand")) {
				strandIdx = i;
			} 
		}
		
		// fail to find fields
		if(Scan.mode.equalsIgnoreCase(Constants.MODE_TARGET)) {
			if(inputSeqIdx == -1 || genomicLociIdx == -1 || strandIdx == -1) {
				System.out.println("Fail to find column names: sequence, location, strand");
				System.out.println("Please specify exact column names in your input file.");
				System.exit(1);
			}
			
		} else {
			if(inputSeqIdx == -1) {
				System.out.println("Fail to find column names: sequence");
				System.out.println("Please specify exact column names in your input file.");
				System.exit(1);
			}
		}
			
		
		if(Scan.mode.equalsIgnoreCase(Constants.MODE_TARGET)) {

			while((line = BR.readLine()) != null) {
				String[] fields = line.split("\t");
				String sequence = fields[inputSeqIdx];
				String genomicLoci = fields[genomicLociIdx];
				String strand = fields[strandIdx];
				
				// if the sequence is invalid
				if(!Validation.checkValidSequence(sequence)) {
					continue;
				}
				
				SequenceRecord record = new SequenceRecord();
				record.sequence = sequence;
				record.strand = strand;
				record.location = genomicLoci;
				
				if(Scan.sequence.equalsIgnoreCase(Constants.SEQUENCE_NUCLEOTIDE)) {
					if(strand.charAt(0) == '-') {
						// rc sequence of nucleotide
						record.sequence = Translator.getReverseComplement(record.sequence);
					}
				}
				
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
		} else {
			while((line = BR.readLine()) != null) {
				String[] fields = line.split("\t");
				String sequence = fields[inputSeqIdx];
				
				// if the sequence is invalid
				if(!Validation.checkValidSequence(sequence)) {
					continue;
				}
				
				SequenceRecord record = new SequenceRecord();
				record.sequence = Scan.isILEqual ? sequence.replace("I", "L") : sequence;
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
		
		// set longestLengthOfInputSequences
		for(SequenceRecord record : records) {
			Scan.longestSequenceLen = Math.max(record.sequence.length(), Scan.longestSequenceLen);
		}
		System.out.println("Records without duplication: "+records.size());
		System.out.println("Longest length of input sequences: "+Scan.longestSequenceLen);
		
		
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
	public static void writeMainOutput (ArrayList<SequenceRecord> records, String baseOutputPath, LocTable locTable) throws IOException {
		writeLibSize(new File(baseOutputPath+".libsize.tsv"));
		BufferedWriter BW = new BufferedWriter(new FileWriter(baseOutputPath+"."+Scan.mode+".tsv"));
		BufferedWriter BWMiss = new BufferedWriter(new FileWriter(baseOutputPath+".miss.tsv"));
		
		// write header
		BW.append(SequenceRecord.header +
				"\t" + Constants.MATCHED_LOCATION +
				"\t" + Constants.MATCHED_MUTATIONS +
				"\t" + Constants.MATCHED_STRAND +
				"\t" + Constants.MATCHED_PEPTIDE +
				"\t" + Constants.MATCHED_NUCLEOTIDE +
				"\t" + Constants.MATCHED_REFNUCLEOTIDE);
		if(Scan.isSingleCellMode) {
			// append barcode ids in whitelist
			for(String barcodeId : BarcodeTable.barcodeIds) {
				BW.append("\t").append(barcodeId);
			}
		} else {
			BW.append("\t" + Constants.MATCHED_READ_COUNT +
					"\t" + Constants.MATCHED_RPHM);
		}
		BW.newLine();
		
		BWMiss.append(SequenceRecord.header);
		BWMiss.newLine();
		
		for(int i=0; i<records.size(); i++) {
			SequenceRecord record = records.get(i);
			String sequence = record.sequence;
			
			ArrayList<LocationInformation> locations = locTable.getLocations(sequence);
			
			for(int j=0; j<record.records.size(); j++) {
				if(locations.size() == 0) {
					BWMiss.append(record.records.get(j));
					BWMiss.newLine();
				} else {
					for(LocationInformation location : locations) {
						// full information (including genomic sequence)
						if(Scan.mode.equalsIgnoreCase(Constants.MODE_TARGET) && 
							!record.location.equalsIgnoreCase(location.location)) {
							continue;
						}
						BW.append(record.records.get(j)).append("\t"+location.getRes());
						BW.newLine();
					}
				}
			}
		}
		BW.close();
		BWMiss.close();
	}
	
	public static void writeLocationLevelOutput (ArrayList<SequenceRecord> records, String baseOutputPath, LocTable locTable) throws IOException {
		BufferedWriter BW = new BufferedWriter(new FileWriter(baseOutputPath+".gloc.tsv"));
		
		// define field index
		BW.append(Constants.MATCHED_PEPTIDE +
				"\t" + Constants.MATCHED_LOCATION +
				"\t" + Constants.MATCHED_STRAND);
		if(Scan.isSingleCellMode) {
			// append barcode ids in whitelist
			for(String barcodeId : BarcodeTable.barcodeIds) {
				BW.append("\t").append(barcodeId);
			}
		} else {
			BW.append("\t" + Constants.MATCHED_READ_COUNT +
					"\t" + Constants.MATCHED_RPHM);
		}
		BW.newLine();
		// end of header //
		
		// write records
		// unique observed sequence.
		Hashtable<String, Hashtable<String, Long>> readCountsTupleLevel = new Hashtable<String, Hashtable<String, Long>>();
		
		for(int i=0; i<records.size(); i++) {
			SequenceRecord record = records.get(i);
			String sequence = record.sequence;
			ArrayList<LocationInformation> locations = locTable.getLocations(sequence);
			
			for(LocationInformation location : locations) {
				if(Scan.mode.equalsIgnoreCase(Constants.MODE_TARGET) && 
						!record.location.equalsIgnoreCase(location.location)) {
						continue;
					}
				
				Hashtable<String, Long> readCounts = location.readCounts;
				// it must be calculated once!
				// peptide level count
				
				String tupleKey = location.obsPeptide+"\t"+location.location+"\t"+location.strand;
				Hashtable<String, Long> sumReads = readCountsTupleLevel.get(tupleKey);
				if(sumReads == null) {
					sumReads = new Hashtable<String, Long>();
				}
				
				Iterator<String> keys = (Iterator<String>) readCounts.keys();
				while(keys.hasNext()) {
					String barcodeId = keys.next();
					Long val = sumReads.get(barcodeId);
					if(val == null) {
						val = 0L;
					}
					sumReads.put(barcodeId, (readCounts.get(barcodeId) + val));
				}
				
				readCountsTupleLevel.put(tupleKey, sumReads);
			}
			
		}
		
		readCountsTupleLevel.forEach((tupleKey, reads)->{
			try {
				if(Scan.isSingleCellMode) {
					BW.append(tupleKey);
					for(String barcodeId : BarcodeTable.barcodeIds) {
						Long read = reads.get(barcodeId);
						if(read == null) {
							read = 0L;
						}
						BW.append("\t"+read);
					}
				} else {
					Long read = reads.get(Constants.DEFAULT_BARCODE_ID);
					BW.append(tupleKey+"\t"+read+"\t"+Utils.getRPHM((double)read));
				}
				
				BW.newLine();
			}catch(IOException ioe) {
				
			}
 		});
		BW.close();
	}
	
	public static void writePeptideLevelOutput (ArrayList<SequenceRecord> records, String baseOutputPath, LocTable locTable) throws IOException {
		BufferedWriter BW = new BufferedWriter(new FileWriter(baseOutputPath+".peptide.tsv"));
		
		// define field index
		BW.append(Constants.MATCHED_PEPTIDE + "(" + Scan.union + ")" +
				"\t" + Constants.MATCHED_NUM_LOCATION);
		if(Scan.isSingleCellMode) {
			// append barcode ids in whitelist
			for(String barcodeId : BarcodeTable.barcodeIds) {
				BW.append("\t").append(barcodeId);
			}
		} else {
			BW.append("\t" + Constants.MATCHED_READ_COUNT +
					"\t" + Constants.MATCHED_RPHM);
		}
		BW.newLine();
		// end of header //
		
		// write records
		// unique observed sequence.
		Hashtable<String, Hashtable<String, Long>> readCountsPeptLevel = new Hashtable<String, Hashtable<String, Long>>();
		Hashtable<String, Hashtable<String, Boolean>> locationsPeptLevel = new Hashtable<String, Hashtable<String, Boolean>>();
		
		for(int i=0; i<records.size(); i++) {
			SequenceRecord record = records.get(i);
			String sequence = record.sequence;
			ArrayList<LocationInformation> locations = locTable.getLocations(sequence);
			
			for(LocationInformation location : locations) {
				if(Scan.mode.equalsIgnoreCase(Constants.MODE_TARGET) && 
						!record.location.equalsIgnoreCase(location.location)) {
						continue;
					}
				
				Hashtable<String, Long> readCounts = location.readCounts;
				// it must be calculated once!
				// peptide level count
				
				Hashtable<String, Long> unionReads = readCountsPeptLevel.get(location.obsPeptide);
				if(unionReads == null) {
					unionReads = new Hashtable<String, Long>();
				}
				
				Iterator<String> keys = (Iterator<String>) readCounts.keys();
				while(keys.hasNext()) {
					String barcodeId = keys.next();
					Long val = unionReads.get(barcodeId);
					if(val == null) {
						val = 0L;
					}
					if(Scan.union.equalsIgnoreCase(Constants.UNION_MAX)) {
						unionReads.put(barcodeId, Math.max(readCounts.get(barcodeId), val));
					} else if(Scan.union.equalsIgnoreCase(Constants.UNION_SUM)){
						unionReads.put(barcodeId, (readCounts.get(barcodeId) + val));
					}
				}
				
				readCountsPeptLevel.put(location.obsPeptide, unionReads);
				
				Hashtable<String, Boolean> gLocationMap = locationsPeptLevel.get(location.obsPeptide);
				if(gLocationMap == null) {
					gLocationMap = new Hashtable<String, Boolean>();
					locationsPeptLevel.put(location.obsPeptide, gLocationMap);
				}
				
				gLocationMap.put(location.location, true);
			}
			
		}

		readCountsPeptLevel.forEach((sequence, reads)->{
			try {
				if(Scan.isSingleCellMode) {
					BW.append(sequence+"\t"+locationsPeptLevel.get(sequence));
					for(String barcodeId : BarcodeTable.barcodeIds) {
						Long read = reads.get(barcodeId);
						if(read == null) {
							read = 0L;
						}
						BW.append("\t"+read);
					}
				} else {
					Long read = reads.get(Constants.DEFAULT_BARCODE_ID);
					BW.append(sequence+"\t"+locationsPeptLevel.get(sequence).size()+"\t"+read+"\t"+Utils.getRPHM((double)read));
				}
				BW.newLine();
			}catch(IOException ioe) {
				
			}
 		});
		
		BW.close();
		
		WriteStatistics.write(baseOutputPath+".stat.tsv", records, readCountsPeptLevel);
	}
	
	private static void writeLibSize (File file) throws IOException {
		BufferedWriter BW = new BufferedWriter(new FileWriter(file));
		BW.append(Scan.libSize+"");
		BW.newLine();
		BW.close();
	}
}
