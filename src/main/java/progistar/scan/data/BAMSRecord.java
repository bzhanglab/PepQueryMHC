package progistar.scan.data;

import java.util.ArrayList;
import java.util.Hashtable;

import org.ahocorasick.trie.Trie;

import progistar.scan.run.Scan;

public class BAMSRecord {

	public static String header = null;
	public static String fileName = null;
	
	public ArrayList<String> records = new ArrayList<String>();
	public String sequence;
	public String strand;
	public String location;
	
	public int start = 0;
	public int end = 0;
	public String chr = null;
	
	public int readCnt = 0;
	
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
	public static Trie getTrie (ArrayList<BAMSRecord> records) {
		ArrayList<String> sequences = new ArrayList<String>();
		Hashtable<String, String> rmDups = new Hashtable<String, String>();
		
		
		for(BAMSRecord record : records) {
			String sequence = record.sequence;
			if(Scan.isILEqual) {
				sequence = sequence.replace("I", "L");
			}
			
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
