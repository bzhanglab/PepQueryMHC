package progistar.scan.function;

import progistar.scan.data.LibraryTable;

public class Utils {

	private static double normaliedValueHM = Math.pow(10, 8);
	private static double normaliedValueHT = Math.pow(10, 5);
	
	public static double getRPHM (double read, String barcode) {
		assert LibraryTable.table.get(barcode) != null;
		
		return (read/LibraryTable.table.get(barcode)) * normaliedValueHM;
	}
	
	public static double getRPHT (double read, String barcode) {
		assert LibraryTable.table.get(barcode) != null;
		
		return (read/LibraryTable.table.get(barcode)) * normaliedValueHT;
	}
	
	
}
