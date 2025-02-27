package progistar.scan.fileIO;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Hashtable;
import java.util.Iterator;
import java.util.LinkedList;

import progistar.scan.data.Annotation;
import progistar.scan.data.BarcodeTable;
import progistar.scan.data.Constants;
import progistar.scan.data.LibraryTable;
import progistar.scan.data.LocTable;
import progistar.scan.data.LocationInformation;
import progistar.scan.data.Parameters;
import progistar.scan.data.SequenceRecord;
import progistar.scan.function.Utils;
import progistar.scan.function.WriteStatistics;

public class WriteOutput {

	private WriteOutput() {}
	
	/**
	 * For target mode
	 * 
	 * @param records
	 * @param file
	 * @throws IOException
	 */
	public static void writeMainOutput (ArrayList<SequenceRecord> records, String baseOutputPath, LocTable locTable) throws IOException {
		writeLibSize(new File(baseOutputPath+".libsize.tsv"));
		BufferedWriter BW = new BufferedWriter(new FileWriter(baseOutputPath+"."+Parameters.mode+".tsv"));
		BufferedWriter BWMiss = new BufferedWriter(new FileWriter(baseOutputPath+".miss.tsv"));
		
		// write header
		BW.append(SequenceRecord.header +
				"\t" + Constants.MATCHED_LOCATION +
				"\t" + Constants.MATCHED_MUTATIONS +
				"\t" + Constants.MATCHED_STRAND +
				"\t" + Constants.MATCHED_PEPTIDE +
				"\t" + Constants.MATCHED_NUCLEOTIDE +
				"\t" + Constants.MATCHED_REFNUCLEOTIDE);
		if(Parameters.isSingleCellMode) {
			// append barcode ids in whitelist
			// write a header for raw read count
			for(String barcodeId : BarcodeTable.barcodeIds) {
				BW.append("\t").append("Match_read_").append(barcodeId);
			}
			// write a header for RPHM
			for(String barcodeId : BarcodeTable.barcodeIds) {
				BW.append("\t").append("Match_RPHT_").append(barcodeId);
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
						if(Parameters.mode.equalsIgnoreCase(Constants.MODE_TARGET) && 
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
		if(Parameters.isSingleCellMode) {
			// append barcode ids in whitelist
			// write a header for raw read count
			for(String barcodeId : BarcodeTable.barcodeIds) {
				BW.append("\t").append("Match_read_").append(barcodeId);
			}
			// write a header for RPHM
			for(String barcodeId : BarcodeTable.barcodeIds) {
				BW.append("\t").append("Match_RPHT_").append(barcodeId);
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
				if(Parameters.mode.equalsIgnoreCase(Constants.MODE_TARGET) && 
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
				if(Parameters.isSingleCellMode) {
					BW.append(tupleKey);
					// write raw read counts
					for(String barcodeId : BarcodeTable.barcodeIds) {
						Long read = reads.get(barcodeId);
						if(read == null) {
							read = 0L;
						}
						BW.append("\t"+read);
					}
					// write RPHTs
					for(String barcodeId : BarcodeTable.barcodeIds) {
						Long read = reads.get(barcodeId);
						if(read == null) {
							read = 0L;
						}
						BW.append("\t"+Utils.getRPHT(read, barcodeId));
					}
				} else {
					Long read = reads.get(Constants.DEFAULT_BARCODE_ID);
					BW.append(tupleKey+"\t"+read+"\t"+Utils.getRPHM((double)read, Constants.DEFAULT_BARCODE_ID));
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
		BW.append(Constants.MATCHED_PEPTIDE + "(" + Parameters.union + ")" +
				"\t" + Constants.MATCHED_NUM_LOCATION);
		if(Parameters.isSingleCellMode) {
			// append barcode ids in whitelist
			// write a header for raw read count
			for(String barcodeId : BarcodeTable.barcodeIds) {
				BW.append("\t").append("Match_read_").append(barcodeId);
			}
			// write a header for RPHM
			for(String barcodeId : BarcodeTable.barcodeIds) {
				BW.append("\t").append("Match_RPHT_").append(barcodeId);
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
				if(Parameters.mode.equalsIgnoreCase(Constants.MODE_TARGET) && 
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
					if(Parameters.union.equalsIgnoreCase(Constants.UNION_MAX)) {
						unionReads.put(barcodeId, Math.max(readCounts.get(barcodeId), val));
					} else if(Parameters.union.equalsIgnoreCase(Constants.UNION_SUM)){
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
				if(Parameters.isSingleCellMode) {
					BW.append(sequence+"\t"+locationsPeptLevel.get(sequence).size());
					// write raw read counts
					for(String barcodeId : BarcodeTable.barcodeIds) {
						Long read = reads.get(barcodeId);
						if(read == null) {
							read = 0L;
						}
						BW.append("\t"+read);
					}
					// write RPHTs
					for(String barcodeId : BarcodeTable.barcodeIds) {
						Long read = reads.get(barcodeId);
						if(read == null) {
							read = 0L;
						}
						BW.append("\t"+Utils.getRPHT(read, barcodeId));
					}
				} else {
					Long read = reads.get(Constants.DEFAULT_BARCODE_ID);
					BW.append(sequence+"\t"+locationsPeptLevel.get(sequence).size()+"\t"+read+"\t"+Utils.getRPHM((double)read, Constants.DEFAULT_BARCODE_ID));
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
		BW.append("Barcode\tLibrary_size");
		BW.newLine();
		ArrayList<String> barcodes = new ArrayList<String>();
		if(Parameters.isSingleCellMode) {
			barcodes = BarcodeTable.barcodeIds;
		} else {
			barcodes.add(Constants.DEFAULT_BARCODE_ID);
		}
		
		for(String barcode : barcodes) {
			Double libSize = LibraryTable.table.get(barcode);
			if(libSize == null) {
				libSize = .0;
			}
			BW.append(barcode+"\t"+libSize);
			BW.newLine();
		}
		BW.close();
	}
	
	
	public static void writeAnnotateOutput (ArrayList<SequenceRecord> records, 
			Hashtable<String, LinkedList<Annotation>> allAnnotations) throws IOException {
		
		// write a result
		BufferedWriter BW = new BufferedWriter(new FileWriter(Parameters.outputBaseFilePath+".annotate.tsv"));
		BW.append(SequenceRecord.header)
		.append("\t").append(Constants.ANNOTATION_COUNT)
		.append("\t").append(Constants.GENE_ID)
		.append("\t").append(Constants.GENE_NAME)
		.append("\t").append(Constants.GENE_STRAND)
		.append("\t").append(Constants.GENE_TYPE)
		.append("\t").append(Constants.CLASS_CODE)
		.append("\t").append(Constants.UNIQUE_CLASS_CODE);
		BW.newLine();
		
		
		StringBuilder geneIds				= new StringBuilder();
		StringBuilder geneNames				= new StringBuilder();
		StringBuilder strands				= new StringBuilder();
		StringBuilder geneTypes				= new StringBuilder();
		StringBuilder geneClassCodes		= new StringBuilder();
		StringBuilder uniqueClassCodes		= new StringBuilder();
		
		for(SequenceRecord sRecord : records) {
			LinkedList<Annotation> annotations = allAnnotations.get(sRecord.getKey());
			assert annotations != null;
			
			for(String record : sRecord.records) {
				
				if(Parameters.stretch) {
					for(Annotation annotation : annotations) {
						BW.append(record);
						BW.append("\t"+annotations.size());
						BW.append("\t").append(annotation.getGeneId());
						BW.append("\t").append(annotation.getGeneName());
						BW.append("\t").append(annotation.getStrand());
						BW.append("\t").append(annotation.getGeneType());
						BW.append("\t").append(annotation.getClassCode()); // class_code
						BW.append("\t").append(annotation.getClassCode()); // unique_class_code
						BW.newLine();
					}
				} else {
					Hashtable<String, Boolean> uniqueClassCode = new Hashtable<String, Boolean>();
					
					for(Annotation annotation : annotations) {
						// build gene ids
						if(geneIds.length() > 0) {
							geneIds.append("|");
						}
						geneIds.append(annotation.getGeneId());
						
						// build gene names
						if(geneNames.length() > 0) {
							geneNames.append("|");
						}
						geneNames.append(annotation.getGeneName());
						
						// build strands
						if(strands.length() > 0) {
							strands.append("|");
						}
						strands.append(annotation.getStrand());
						
						// build gene names
						if(geneTypes.length() > 0) {
							geneTypes.append("|");
						}
						geneTypes.append(annotation.getGeneType());
						
						// build gene names
						if(geneClassCodes.length() > 0) {
							geneClassCodes.append("|");
						}
						geneClassCodes.append(annotation.getClassCode());
						
						// build unique class code
						if(uniqueClassCode.get(annotation.getClassCode()) == null) {
							if(uniqueClassCodes.length() > 0) {
								uniqueClassCodes.append("|");
							}
							uniqueClassCodes.append(annotation.getClassCode());
							uniqueClassCode.put(annotation.getClassCode(), true);
						}
					}
					BW.append(record);
					BW.append("\t"+annotations.size());
					BW.append("\t").append(geneIds.toString()); geneIds.setLength(0);
					BW.append("\t").append(geneNames.toString()); geneNames.setLength(0);
					BW.append("\t").append(strands.toString()); strands.setLength(0);
					BW.append("\t").append(geneTypes.toString()); geneTypes.setLength(0);
					BW.append("\t").append(geneClassCodes.toString()); geneClassCodes.setLength(0);
					BW.append("\t").append(uniqueClassCodes.toString()); uniqueClassCodes.setLength(0);
					BW.newLine();
				}
				
				
			}
		}
		
		BW.close();
	}
}
