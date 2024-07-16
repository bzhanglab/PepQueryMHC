package progistar.scan.data;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.Hashtable;

import progistar.scan.run.Scan;

public class BarcodeTable {

	public static Hashtable<String, Boolean> hasBarcode = new Hashtable<String, Boolean>();
	
	public static void load() {
		assert Scan.whitelistFile != null;
		assert Scan.isSingleCellMode;
		
		try {
			BufferedReader BR = new BufferedReader(new FileReader(Scan.whitelistFile));
			String line = null;
			
			BR.readLine(); // skip header
			while((line = BR.readLine()) != null) {
				String[] fields = line.split("\t");
				hasBarcode.put(fields[0], true);
			}
			
			BR.close();
			
			System.out.println("A total of "+hasBarcode.size()+" barcodes were saved.");
		}catch(IOException ioe) {
			System.out.println("Fail to load white list: "+Scan.whitelistFile.getName());
		}
	}
}
