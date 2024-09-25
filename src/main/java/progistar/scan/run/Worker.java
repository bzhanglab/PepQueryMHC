package progistar.scan.run;

import java.util.concurrent.Callable;

import progistar.scan.data.Constants;
import progistar.scan.function.CheckMemory;
import progistar.scan.function.ScanModeRun;
import progistar.scan.function.StrandDetection;
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
		
		// strand auto detection
		if(this.task.type == Constants.TYPE_STRAND_DETECTION_TASK) {
			StrandDetection.runDetection(task);
		} 
		// else: general scan / target modes
		else {
			if(Scan.mode.equalsIgnoreCase(Constants.MODE_SCAN)) {
				ScanModeRun.runScanMode(task);
			} else if(Scan.mode.equalsIgnoreCase(Constants.MODE_TARGET)) {
				TargetModeRun.runTargetMode(task);
			}
		}
		
		
		// call once
		this.task.peakMemory = CheckMemory.checkUsedMemoryMB();
		
		return task.getTaskInfo();
	}
}
