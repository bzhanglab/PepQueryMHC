package progistar.scan.run;

import java.util.concurrent.Callable;

import progistar.scan.data.Constants;
import progistar.scan.function.CheckMemory;
import progistar.scan.function.ScanModeRun;
import progistar.scan.function.TargetModeRun;

public class Worker implements Callable<String> {

	private Task task;
	
	public Worker (Task task) {
		super();
		this.task = task;
	}
	
	public String call() {
		if(Scan.verbose) {
			System.out.println(task.getTaskInfo());
		}
		if(Scan.mode.equalsIgnoreCase(Constants.MODE_SCAN)) {
			ScanModeRun.runScanMode(task);
		} else if(Scan.mode.equalsIgnoreCase(Constants.MODE_TARGET)) {
			TargetModeRun.runTargetMode(task);
		}
		
		// call once
		this.task.peakMemory = CheckMemory.checkUsedMemoryMB();
		
		return task.getTaskInfo();
	}
}
