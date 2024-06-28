package progistar.scan.function;

import progistar.scan.run.Scan;

public class Utils {

	private static double normaliedValue = Math.pow(10, 8);
	public static double getRPHM (double read) {
		return (read/Scan.libSize) * normaliedValue;
	}
}
