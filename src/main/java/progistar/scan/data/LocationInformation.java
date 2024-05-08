package progistar.scan.data;

import java.util.ArrayList;
import java.util.Hashtable;

import progistar.scan.function.Translator;

public class LocationInformation {
	public String inputSequence;
	
	public String location;
	public String mutation;
	public String obsNucleotide;
	public String refNucleotide;
	public String obsPeptide;
	public String refPeptide;
	public long readCount = 1;
	public char strand;
	
	
	public String getKey () {
		return location+"\t"+strand+"\t"+obsNucleotide+"\t"+refNucleotide;
	}
	
	public String getRes () {
		return location+"\t"+mutation+"\t"+strand+"\t"+obsNucleotide+"\t"+obsPeptide+"\t"+refNucleotide+"\t"+readCount;
	}
	
	public void calMetaInfo () {
		this.calMutation();
		assert this.obsNucleotide != null;
		
		if(strand == '+') {
			this.obsPeptide = Translator.translation(this.obsNucleotide, 0);
		} else if(strand == '-') {
			this.obsPeptide = Translator.translation(Translator.getReverseComplement(this.obsNucleotide), 0);
		}
	}
	
	public void calMutation () {
		String[] locations = location.split("\\|");
		int mPos =0;
		ArrayList<Mutation> mutations = new ArrayList<Mutation>();
		for(int i=0; i<locations.length; i++) {
			String location = locations[i];
			
			// unmapped location
			if(location.equalsIgnoreCase("-")) {
				this.mutation = "-";
			} else {
				String chr = location.split("\\:")[0];
				String[] fields = location.split("\\:")[1].split("\\-");
				int start = Integer.parseInt(fields[0]);
				int end = Integer.parseInt(fields[1]);
				
				for(int j=start; j<=end; j++) {
					// mutation
					if(this.obsNucleotide.charAt(mPos) != this.refNucleotide.charAt(mPos)) {
						Mutation mutation = new Mutation();
						mutation.altSeq = this.obsNucleotide.charAt(mPos)+"";
						mutation.chrName = chr;
						mutation.refSeq = this.refNucleotide.charAt(mPos)+"";
						mutation.genomicPosition = j;
						// softclip
						if(this.refNucleotide.charAt(mPos) == '*') {
							mutation.type = Constants.CLP;
						} 
						// insertion
						else if(this.refNucleotide.charAt(mPos) == '-') {
							j--;
							mutation.genomicPosition = j;
							mutation.type = Constants.INS;
						}
						// deletion
						else if(this.obsNucleotide.charAt(mPos) == '-') {
							mutation.type = Constants.DEL;
						}
						// snp
						else {
							mutation.type = Constants.SNP;
						}
						mutations.add(mutation);
					}
					mPos++;
				}
			}
		}
		
		if(mutations.size() == 0) {
			this.mutation = "-";
		} else {
			this.mutation = "";
			for(int i=0; i<mutations.size()-1; i++) {
				Mutation prevM = mutations.get(i);
				Mutation nextM = mutations.get(i+1);
				
				if(prevM.type == Constants.SNP) {
					continue;
				}
				
				if(prevM.type == nextM.type) {
					prevM.refSeq += nextM.refSeq;
					prevM.altSeq += nextM.altSeq;
					mutations.remove(i+1);
					i--;
				}
			}
			
			for(int i=0; i<mutations.size(); i++) {
				// lower to upper
				mutations.get(i).refSeq = mutations.get(i).refSeq.toUpperCase();
				if(mutation.length() != 0) {
					mutation += "|";
				}
				mutation += mutations.get(i).toString();
					
			}
		}
	}
}
