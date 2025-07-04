package progistar.scan.function;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Hashtable;
import java.util.Iterator;

import progistar.scan.data.BarcodeTable;
import progistar.scan.data.Constants;
import progistar.scan.data.Parameters;
import progistar.scan.data.SequenceRecord;

public class WriteStatistics {

	/**
	 * count peptides and calculate proportion of matches.
	 * 
	 * 
	 * @param fileName
	 * @param readCountsPeptLevel
	 * @param readCountsRandomPeptLevel
	 * @throws IOException
	 */
	public static void write (String fileName,
							  ArrayList<SequenceRecord> records,
							  Hashtable<String, Hashtable<String, Long>> readCountsPeptLevel) 
							  throws IOException {
		BufferedWriter BW = new BufferedWriter(new FileWriter(fileName));
		
		long[] wholeTable = new long[Parameters.longestSequenceLen+1];
		long[] matchTable = new long[Parameters.longestSequenceLen+1];
		int maxLen = 0;
		int minLen = Integer.MAX_VALUE;
		
		// count input and decoy sequences
		Hashtable<String, Boolean> checkList1 = new Hashtable<String, Boolean>();
		for(SequenceRecord record : records) {
			String sequence = record.sequence;
			if(Parameters.isILEqual) {
				sequence = sequence.replace("I", "L");
			}
			
			if(checkList1.get(sequence) == null) {
				int sequenceLen = sequence.length();
				wholeTable[sequenceLen]++;
				wholeTable[0]++;
				checkList1.put(sequence, true);
				
				maxLen = Math.max(maxLen, sequenceLen);
				minLen = Math.min(minLen, sequenceLen);
			}
		}
		
		// check list to deal with I/L equal option
		// input and random sequences share the hashtable, assuming that there are no overlaps between them.
		Hashtable<String, Boolean> checkList2 = new Hashtable<String, Boolean>();
		// count for input sequences
		readCountsPeptLevel.forEach((sequence, reads) -> {
			if(Parameters.isILEqual) {
				sequence = sequence.replace("I", "L");
			}
			
			Long sum = 0L;
			if(Parameters.isSingleCellMode) {
				
				// calculate sum of reads
				for(String barcodeId : BarcodeTable.barcodeIds) {
					Long read = reads.get(barcodeId);
					if(read == null) {
						read = 0L;
					}
					sum += read;
				}
			} else {
				sum = reads.get(Constants.DEFAULT_BARCODE_ID);
			}
			
			if(checkList2.get(sequence) == null && sum > 0L) {
				int sequenceLen = sequence.length();
				matchTable[sequenceLen]++;
				matchTable[0]++;
				checkList2.put(sequence, true);
			}
		});
		
		BW.append("Length\tInput_sequences\tMatched_sequences\tProportion");
		BW.newLine();
		
		// for input
		for(int i=0; i<wholeTable.length; i++) {
			if(i == 0 || (i >= minLen && i <= maxLen)) {
				if(i == 0) {
					BW.append("Total\t");
				} else {
					BW.append(i+"\t");
				}
				if(wholeTable[i] > 0) {
					BW.append(wholeTable[i]+"\t")
					.append(matchTable[i]+"\t")
					.append(""+((double)(matchTable[i])) / ((double)(wholeTable[i])));
				} else {
					BW.append(wholeTable[i]+"\t")
					.append(matchTable[i]+"\t")
					.append("NaN");
				}
				BW.newLine();	
			}
			
		}
		
		BW.close();
	}
	
}
