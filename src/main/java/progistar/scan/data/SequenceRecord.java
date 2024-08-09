package progistar.scan.data;

import java.util.ArrayList;
import java.util.Hashtable;

import org.ahocorasick.trie.Trie;

import progistar.scan.run.Scan;

public class SequenceRecord {

	public static String header = null;
	public static String fileName = null;
	
	// save duplicated records (based on key)
	public ArrayList<String> records = new ArrayList<String>();
	public String sequence;
	public String strand;
	public String location;
	
	public int start = 0;
	public int end = 0;
	public String chr = null;
	
	public Hashtable<String, Long> readCounts = new Hashtable<String, Long>();
	
	public String getKey () {
		return (sequence+"_"+location+"_"+strand);
	}
	

	/**
	 * If mode is all, than it uses all records.
	 * On the other hand, it uses only records with unmapped reads.
	 * 
	 * 
	 * @param records
	 * @param mode
	 * @return
	 */
	public static Trie getTrie (ArrayList<SequenceRecord> records) {
		ArrayList<String> sequences = new ArrayList<String>();
		Hashtable<String, String> rmDups = new Hashtable<String, String>();
		
		for(SequenceRecord record : records) {
			String sequence = record.sequence;
			if(rmDups.get(sequence) == null) {
				sequences.add(sequence);
				rmDups.put(sequence, "");
			}
		}
		
		if(sequences.size() == 0) {
			return null;
		} else {
			return Trie.builder().addKeywords(sequences).build();
		}
	}
}
