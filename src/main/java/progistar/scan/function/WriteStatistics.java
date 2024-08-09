package progistar.scan.function;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Hashtable;

import progistar.scan.data.SequenceRecord;
import progistar.scan.run.Scan;

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
							  Hashtable<String, Hashtable<String, Long>> readCountsPeptLevel, 
							  Hashtable<String, Hashtable<String, Long>> readCountsRandomPeptLevel) 
							  throws IOException {
		BufferedWriter BW = new BufferedWriter(new FileWriter(fileName));
		
		long[][] wholeTable = new long[Scan.longestSequenceLen+1][2];
		long[][] matchTable = new long[Scan.longestSequenceLen+1][2];
		
		// count input and decoy sequences
		for(SequenceRecord record : records) {
			String sequence = record.sequence;
			int length = sequence.length();
			
			if(record.isRandom) {
				wholeTable[length][1]++;
				wholeTable[0][1]++;
			} else {
				wholeTable[length][0]++;
				wholeTable[0][0]++;
			}
		}
		
		// check list to deal with I/L equal option
		// input and random sequences share the hashtable, assuming that there are no overlaps between them.
		Hashtable<String, Boolean> checkList = new Hashtable<String, Boolean>();
		// count for input sequences
		readCountsPeptLevel.forEach((sequence, cnt) -> {
			if(Scan.isILEqual) {
				sequence = sequence.replace("I", "L");
			}
			
			if(checkList.get(sequence) == null) {
				int sequenceLen = sequence.length();
				matchTable[sequenceLen][0]++;
				checkList.put(sequence, true);
			}
		});
		// count for input decoy sequences
		readCountsRandomPeptLevel.forEach((sequence, cnt) -> {
			if(Scan.isILEqual) {
				sequence = sequence.replace("I", "L");
			}
			
			if(checkList.get(sequence) == null) {
				int sequenceLen = sequence.length();
				matchTable[sequenceLen][1]++;
				checkList.put(sequence, true);
			}
		});
		
		for(int i=0; i<matchTable.length; i++) {
			matchTable[0][0] += matchTable[i][0];
			matchTable[0][1] += matchTable[i][1];
		}
		
		BW.append("Label\tLength\tNumOfSequences\tNumOfMatches\tProportion");
		BW.newLine();
		
		// for input
		for(int i=0; i<wholeTable.length; i++) {
			if(wholeTable[i][0] > 0) {
				BW.append("Input\t")
				.append(i+"\t")
				.append(wholeTable[i][0]+"\t")
				.append(matchTable[i][0]+"\t")
				.append(""+((double)(matchTable[i][0])) / ((double)(wholeTable[i][0])));
				BW.newLine();
			}
		}
		
		// for random
		for(int i=0; i<wholeTable.length; i++) {
			if(wholeTable[i][1] > 0) {
				BW.append("Random\t")
				.append(i+"\t")
				.append(wholeTable[i][1]+"\t")
				.append(matchTable[i][1]+"\t")
				.append(""+((double)(matchTable[i][1])) / ((double)(wholeTable[i][1])));
				BW.newLine();
			}
		}
		
		
		BW.close();
	}
	
}
