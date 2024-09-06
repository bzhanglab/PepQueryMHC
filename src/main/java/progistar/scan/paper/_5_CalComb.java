package progistar.scan.paper;

import java.io.IOException;
import java.util.Hashtable;

public class _5_CalComb {

	public static Hashtable<Character, Double> loadTable(boolean isILEqual) {
		Hashtable<Character, Double> table = new Hashtable<Character, Double>();
		
		table.put('M', 1.);
		table.put('W', 1.);
		
		table.put('F', 2.);
		table.put('Y', 2.);
		table.put('H', 2.);
		table.put('Q', 2.);
		table.put('N', 2.);
		table.put('K', 2.);
		table.put('D', 2.);
		table.put('E', 2.);
		table.put('C', 2.);
		
		table.put('I', 3.);
		
		table.put('V', 4.);
		table.put('P', 4.);
		table.put('T', 4.);
		table.put('A', 4.);
		table.put('G', 4.);
		
		table.put('L', 6.);
		table.put('S', 6.);
		table.put('R', 6.);
		
		// for I = L
		if(isILEqual) {
			double sum = table.get('I') +table.get('L');
			table.put('L', sum);
			table.put('I', sum);
		}
		
		
		return table;
	}
	
	public static double getCombi (String peptide, Hashtable<Character, Double> table) {
		double cnt = 1;
		
		for(int i=0; i<peptide.length(); i++) {
			char aa = peptide.charAt(i);
			cnt *= table.get(aa);
		}
		
		return cnt;
	}
	
	public static void main(String[] args) throws IOException {
		Hashtable<Character, Double> table = loadTable(false);
		// check
		double[] sum = new double[1];
		table.forEach((AA, val)->{
			sum[0] += val;
		});
		System.out.println(sum[0]);
		
		double total = 0;
		String[] list = {"LPLQAQPSA", "APLRAGWAA", "TTLLQGSLVVK",
				"EVYHTTVLK", "SGFSFQVTIRK", "TALLWSLRK", "TTALLWSLRK", "SSSDLKYLK", 
				"FSFQVTIRK", "AVVKAEQHLLK"};
		for(String l : list) {
			double c = getCombi(l, table);
			System.out.println(l+"\t"+c);
			total += c;
		}
		System.out.println("Total\t"+total);
		
	}
}
