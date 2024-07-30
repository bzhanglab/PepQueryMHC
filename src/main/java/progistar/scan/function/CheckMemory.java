package progistar.scan.function;

public class CheckMemory {

	public static long checkUsedMemoryByte () {
		long totalMemory = Runtime.getRuntime().totalMemory();
		long freeMemory = Runtime.getRuntime().freeMemory();
		
		return totalMemory - freeMemory;
	}
	
	public static long checkUsedMemoryMB () {
		long usedMemoryByte = checkUsedMemoryByte();
		
		return usedMemoryByte / (1024*1024);
	}
}
