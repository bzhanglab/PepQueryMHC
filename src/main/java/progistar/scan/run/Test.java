package progistar.scan.run;

import java.io.IOException;

import progistar.scan.data.Phred;

public class Test {

	public static void main(String[] args) throws IOException {
		Phred.loadTable();
		
		String target1 = "JJ<JJJJJJJJJJJFJA7A-AJJAAJJ";
		String target2 = "F-<7AJJJ77-AA<F<<7F<AFAJFJJ";
		
		System.out.println(Phred.getProbOfAtLeastOneError(target1));
		System.out.println(Phred.getProbOfAtLeastOneError(target2));
		
		String[] a = new String[3];
		
		System.out.println(a[0]+"\t"+a[1]);
	}
}
