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
				String sequence = fields[obsSeqIdx];
				
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
			
			if(Scan.isRandom) {
				indexedRecords.forEach((key, record)->{
					SequenceRecord rRecord = new SequenceRecord();
					rRecord.sequence = Random.getReverseSequence(record.sequence);
					rRecord.strand = Constants.NULL;
					rRecord.location = Constants.NULL;
					rRecord.isRandom = true;
					
					if(indexedRecords.get(rRecord.getKey()) == null) {
						records.add(rRecord);
					}
				});
				
				int numOfRandomSequences = 0;
				for(SequenceRecord record : records) {
					if(record.isRandom) {
						numOfRandomSequences ++;
					}
				}
				System.out.println("The number of "+numOfRandomSequences+" random sequences were generated.");
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
	public static void writeRecords (ArrayList<SequenceRecord> records, File file) throws IOException {
		writeLibSize(file);
		
		BufferedWriter BW = new BufferedWriter(new FileWriter(file));
		BufferedWriter BWPeptCount = new BufferedWriter(new FileWriter(file.getAbsolutePath()+".pept_count.tsv"));
		
		// write header
		if(Scan.isSingleCellMode) {
			BW.append(SequenceRecord.header);
			BWPeptCount.append("Sequence");
			
			// append barcode ids in whitelist
			for(String barcodeId : BarcodeTable.barcodeIds) {
				BW.append("\t").append(barcodeId);
				BWPeptCount.append("\t").append(barcodeId);
			}
			
		} else {
			BW.append(SequenceRecord.header+"\tReadCount\tRPHM");
			BWPeptCount.append("Sequence\tReadCount\tRPHM");
		}
		
		BW.newLine();
		BWPeptCount.newLine();
		
		Hashtable<String, Hashtable<String, Long>> readCountsPeptLevel = new Hashtable<String, Hashtable<String, Long>>();
		
		// write records
		for(int i=0; i<records.size(); i++) {
			SequenceRecord record = records.get(i);
			
			Hashtable<String, Long> readCounts = record.readCounts;
			for(int j=0; j<record.records.size(); j++) {
				if(Scan.isSingleCellMode) {
					BW.append(record.records.get(j));
					for(String barcodeId : BarcodeTable.barcodeIds) {
						Long read = readCounts.get(barcodeId);
						if(read == null) {
							read = 0L;
						}
						BW.append("\t"+read);
					}
				} else {
					Long read = readCounts.get(Constants.DEFAULT_BARCODE_ID);
					if(read == null) {
						read = 0L;
					}
					BW.append(record.records.get(j)).append("\t"+read+"\t"+Utils.getRPHM((double)read));
				}
				BW.newLine();
			}
			
			Hashtable<String, Long> sumReads = readCountsPeptLevel.get(record.sequence);
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
			
			readCountsPeptLevel.put(record.sequence, sumReads);
		}
		
		readCountsPeptLevel.forEach((sequence, readCounts)->{
			try {
				if(Scan.isSingleCellMode) {
					BWPeptCount.append(sequence);
					for(String barcodeId : BarcodeTable.barcodeIds) {
						Long read = readCounts.get(barcodeId);
						if(read == null) {
							read = 0L;
						}
						BWPeptCount.append("\t"+read);
					}
				} else {
					Long read = readCounts.get(Constants.DEFAULT_BARCODE_ID);
					BWPeptCount.append(sequence+"\t"+read+"\t"+Utils.getRPHM((double)read));
				}
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
	public static void writeRecords (ArrayList<SequenceRecord> records, File file, ArrayList<Task> tasks) throws IOException {
		writeLibSize(file);
		
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
		if(Scan.isSingleCellMode) {
			BW.append(SequenceRecord.header+"\tLocation\tMutations\tStrand\tObsNucleotide\tObsPeptide\tRefNucleotide");
			BWGenomicTuple.append("ObsPeptide\tLocation\tStrand");
			BWPeptCount.append("ObsPeptide\tNumLocations");
			
			// append barcode ids in whitelist
			for(String barcodeId : BarcodeTable.barcodeIds) {
				BW.append("\t").append(barcodeId);
				BWGenomicTuple.append("\t").append(barcodeId);
				BWPeptCount.append("\t").append(barcodeId);
			}
			
		} else {
			BW.append(SequenceRecord.header+"\tLocation\tMutations\tStrand\tObsNucleotide\tObsPeptide\tRefNucleotide\tReadCount\tRPHM");
			BWGenomicTuple.append("ObsPeptide\tLocation\tStrand\tReadCount\tRPHM");
			BWPeptCount.append("ObsPeptide\tNumLocations\tReadCount\tRPHM");
		}
		BW.newLine();
		BWGenomicTuple.newLine();
		
		BWNotFound.append(SequenceRecord.header+"\tLocation");
		BWNotFound.newLine();
		BWPeptCount.newLine();
		
		
		// write records
		// unique observed sequence.
		Hashtable<String, Hashtable<String, Long>> readCountsPeptLevel = new Hashtable<String, Hashtable<String, Long>>();
		Hashtable<String, Integer> locationsPeptLevel = new Hashtable<String, Integer>();
		
		Hashtable<String, Hashtable<String, Long>> readCountsRandomPeptLevel = new Hashtable<String, Hashtable<String, Long>>();
		Hashtable<String, Integer> locationsRandomPeptLevel = new Hashtable<String, Integer>();
		
		Hashtable<String, Hashtable<String, Long>> readCountsTupleLevel = new Hashtable<String, Hashtable<String, Long>>();
		Hashtable<String, Boolean> isUniqueCal = new Hashtable<String, Boolean>();
		
		for(int i=0; i<records.size(); i++) {
			SequenceRecord record = records.get(i);
			String sequence = record.sequence;
			ArrayList<LocationInformation> locations = locTable.getLocations(sequence);
			
			// the records must be an unique item!
			if(isUniqueCal.get(sequence) != null) {
				System.err.println("Severe: the records is not unique!");
			}
			isUniqueCal.put(sequence, true);
			
			for(LocationInformation location : locations) {
				Hashtable<String, Long> readCounts = location.readCounts;
				// it must be calculated once!
				// peptide level count
				
				// random count
				if(record.isRandom) {
					Hashtable<String, Long> sumReads = readCountsRandomPeptLevel.get(location.obsPeptide);
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
					
					readCountsRandomPeptLevel.put(location.obsPeptide, sumReads);
					locationsRandomPeptLevel.put(location.obsPeptide, locations.size());
				} 
				// non-random count
				else {
					Hashtable<String, Long> sumReads = readCountsPeptLevel.get(location.obsPeptide);
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
					
					readCountsPeptLevel.put(location.obsPeptide, sumReads);
					locationsPeptLevel.put(location.obsPeptide, locations.size());
					
					// tuple level count
					String tupleKey = location.obsPeptide+"\t"+location.location+"\t"+location.strand;
					sumReads = readCountsTupleLevel.get(tupleKey);
					if(sumReads == null) {
						sumReads = new Hashtable<String, Long>();
					}
					
					keys = (Iterator<String>) readCounts.keys();
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
			
			
			
			// if there are duplicated records, then this size > 1
			// if there is no duplication, then this size = 1
			
			// pass random sequence
			// random sequences are only written in the peptide count.
			if(!record.isRandom) {
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
			
		}
		BWNotFound.close();
		BW.close();
		
		readCountsPeptLevel.forEach((sequence, reads)->{
			try {
				if(Scan.isSingleCellMode) {
					BWPeptCount.append(sequence+"\t"+locationsPeptLevel.get(sequence));
					for(String barcodeId : BarcodeTable.barcodeIds) {
						Long read = reads.get(barcodeId);
						if(read == null) {
							read = 0L;
						}
						BWPeptCount.append("\t"+read);
					}
				} else {
					Long read = reads.get(Constants.DEFAULT_BARCODE_ID);
					BWPeptCount.append(sequence+"\t"+locationsPeptLevel.get(sequence)+"\t"+read+"\t"+Utils.getRPHM((double)read));
				}
				BWPeptCount.newLine();
			}catch(IOException ioe) {
				
			}
 		});
		
		BWPeptCount.close();
		
		readCountsTupleLevel.forEach((tupleKey, reads)->{
			try {
				if(Scan.isSingleCellMode) {
					BWGenomicTuple.append(tupleKey);
					for(String barcodeId : BarcodeTable.barcodeIds) {
						Long read = reads.get(barcodeId);
						if(read == null) {
							read = 0L;
						}
						BWGenomicTuple.append("\t"+read);
					}
				} else {
					Long read = reads.get(Constants.DEFAULT_BARCODE_ID);
					BWGenomicTuple.append(tupleKey+"\t"+read+"\t"+Utils.getRPHM((double)read));
				}
				
				BWGenomicTuple.newLine();
			}catch(IOException ioe) {
				
			}
 		});
		BWGenomicTuple.close();
		
		
		
		// if calculate random distribution is on :
		if(Scan.isRandom) {
			BufferedWriter BWRandomPeptCount = new BufferedWriter(new FileWriter(file.getAbsolutePath()+".pept_count.random.tsv"));
			if(Scan.isSingleCellMode) {
				BWRandomPeptCount.append("rObsPeptide\tNumLocations");
				for(String barcodeId : BarcodeTable.barcodeIds) {
					BWRandomPeptCount.append("\t").append(barcodeId);
				}
			} else {
				BWRandomPeptCount.append("rObsPeptide\tNumLocations\tReadCount\tRPHM");
			}
			BWRandomPeptCount.newLine();
			readCountsRandomPeptLevel.forEach((sequence, reads)->{
				try {
					if(Scan.isSingleCellMode) {
						BWRandomPeptCount.append(sequence+"\t"+locationsRandomPeptLevel.get(sequence));
						for(String barcodeId : BarcodeTable.barcodeIds) {
							Long read = reads.get(barcodeId);
							if(read == null) {
								read = 0L;
							}
							BWRandomPeptCount.append("\t"+read);
						}
					} else {
						Long read = reads.get(Constants.DEFAULT_BARCODE_ID);
						BWRandomPeptCount.append(sequence+"\t"+locationsRandomPeptLevel.get(sequence)+"\t"+read+"\t"+Utils.getRPHM((double)read));
					}
					BWRandomPeptCount.newLine();
				}catch(IOException ioe) {
					
				}
	 		});
			
			BWRandomPeptCount.close();
			
		}
	}
	
	private static void writeLibSize (File file) throws IOException {
		BufferedWriter BW = new BufferedWriter(new FileWriter(file.getAbsolutePath()+".libsize"));
		BW.append(Scan.libSize+"");
		BW.close();
	}
}
