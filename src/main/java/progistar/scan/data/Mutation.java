package progistar.scan.data;

public class Mutation {

	public byte type;
	
	public String altSeq;
	public String refSeq;
	public int relPos;
	
	public int genomicPosition;
	public String chrName;
	
	public String toString () {
		if(type == Constants.SNP) return chrName +":" +genomicPosition+refSeq+">"+altSeq;
		else if(type == Constants.INS){
			return chrName +":" +genomicPosition+"ins"+altSeq;
		} else if(type == Constants.DEL){
			return chrName +":" +genomicPosition+"del"+refSeq;
		} else if(type == Constants.CLP) {
			return chrName +":" +genomicPosition+"clp"+refSeq;
		}
		
		return "NA";
	}
}
