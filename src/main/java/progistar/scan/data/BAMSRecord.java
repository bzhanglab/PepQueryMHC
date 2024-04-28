package progistar.scan.data;

import java.util.ArrayList;
import java.util.Hashtable;

import org.ahocorasick.trie.Trie;

public class BAMSRecord {

	public static String header = null;
	
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
			if(rmDups.get(record.sequence) == null) {
				sequences.add(record.sequence);
				rmDups.put(record.sequence, "");
			}
		}
		
		if(sequences.size() == 0) {
			return null;
		} else {
			return Trie.builder().addKeywords(sequences).build();
		}
	}
}
