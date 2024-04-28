package progistar.scan.run;

import java.util.concurrent.Callable;

import progistar.scan.data.Constants;
import progistar.scan.function.ScanModeRun;
import progistar.scan.function.TargetModeRun;

public class Worker implements Callable<String> {

	private Task task;
	private int workerID = -1;
	
	public Worker (int workerID, Task task) {
		super();
		this.task = task;
		this.workerID = workerID;
	}
	
	public String call() {
		System.out.println(task.getTaskInfo()+" by "+this.workerID);
		if(Scan.mode.equalsIgnoreCase(Constants.MODE_SCAN)) {
			ScanModeRun.runScanMode(task);
		} else if(Scan.mode.equalsIgnoreCase(Constants.MODE_TARGET)) {
			TargetModeRun.runTargetMode(task);
		}
		
		return task.getTaskInfo();
	}
}
