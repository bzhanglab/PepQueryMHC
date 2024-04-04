package progistar.scan.run;

import progistar.scan.data.Constants;
import progistar.scan.function.TargetModeRun;

public class Worker extends Thread {

	private Task task;
	private int workerID = -1;
	
	public Worker (int workerID, Task task) {
		super();
		this.task = task;
		this.workerID = workerID;
	}
	
	public void run() {
		System.out.println(task.getTaskInfo()+" by "+this.workerID);
		if(Scan.mode.equalsIgnoreCase(Constants.MODE_FULL)) {
			
		} else if(Scan.mode.equalsIgnoreCase(Constants.MODE_TARGET)) {
			TargetModeRun.runTargetMode(task);
		}
	}
}
