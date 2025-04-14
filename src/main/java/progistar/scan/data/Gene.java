package progistar.scan.data;

import java.util.Hashtable;
import java.util.Iterator;

import progistar.scan.function.IndexConvertor;

public class Gene implements Comparable<Gene> {

	public int chrIdx;
	public int start;
	public int end;
	
	
	public int min = -1;
	public int max = -1;
	
	public String id;
	public String name;
	public String type;
	
	public Hashtable<String, Transcript> transcripts = new Hashtable<String, Transcript>();
	
	// this value represents a visit status of gene.
	// it is automatically increased by "GeneArray"
	// do not handle this value directly.
	public int mark = -1;
	
	@Override
	public int compareTo(Gene o) {
		if(this.start < o.start) {
			return -1;
		} else if(this.start > o.start) {
			return 1;
		}
		return 0;
	}
	
	public String toString() {
		return id+" ("+name+"), "+IndexConvertor.indexToChr(this.chrIdx)+", "+this.type+" :"+start+"-"+end;
	}
	
	public Annotation[] annotate (SequenceRecord sRecord) {
		int size = transcripts.size();
		Annotation[] annotates = new Annotation[size];
		
		Iterator<String> iterators = (Iterator<String>) transcripts.keys();
		
		int idx = 0;
		while(iterators.hasNext()) {
			String id = iterators.next();
			annotates[idx] = transcripts.get(id).getClassCode(sRecord);
			annotates[idx].gene = this;
			idx++;
		}
		
		return annotates;
	}
}
