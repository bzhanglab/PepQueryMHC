package progistar.scan.data;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.Hashtable;

public class LibraryTable {

	public static Hashtable<String, Double> table = new Hashtable<String, Double>();
	
	public static void loadTable (File file) throws IOException {
		System.out.println("Load library table...");
		BufferedReader BR = new BufferedReader(new FileReader(file));
		String line = null;
		
		// the header should contain "Barcode" and "Library_size"
		BR.readLine(); // skip header
		int total = 0;
		while((line = BR.readLine()) != null) {
			String[] fields = line.split("\t");
			String barcode = fields[0];
			Double libSize = Double.parseDouble(fields[1]);
			
			table.put(barcode, libSize);
			total ++;
		}
		
		if(table.size() != total) {
			System.out.println("There were duplicated barcodes: "+(total - table.size()));
			System.out.println("This may lead to a wrong RPHM calculation...!");
		}
		
		BR.close();
	}
	
	public static boolean isEmpty() {
		if(table.size() == 0) {
			return true;
		}
		
		return false;
	}
}
