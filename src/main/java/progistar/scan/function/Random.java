package progistar.scan.function;

import java.util.ArrayList;

public class Random {

	
	public static String getReverseSequence (String str) {
		return new StringBuilder(str).reverse().toString();
	}
	
	
	public static ArrayList<String> getRandomSequences (String str) {
		ArrayList<String> list = new ArrayList<String>();
		String revSequence = getReverseSequence(str);
		list.add(revSequence);
		
		int len = str.length();
		// forward shifted
		//for(int i=1; i<len; i++) {
		//	list.add(str.substring(len-i, len) + str.substring(0, len-i));
		//}
		// reverse shifted
		for(int i=1; i<len; i++) {
			list.add(revSequence.substring(len-i, len) + revSequence.substring(0, len-i));
		}
		
		return list;
	}
}
