package progistar.scan.run;

import java.util.concurrent.Callable;

import progistar.scan.data.Constants;
import progistar.scan.data.Parameters;
import progistar.scan.function.CheckMemory;
import progistar.scan.function.ScanModeRun;
import progistar.scan.function.StrandDetection;
import progistar.scan.function.TargetModeRun;

public class Worker implements Callable<String> {

	private static int done = 0;
	private static int totalTasks = 0;
	private Task task;
	
	/**
	 * Use it when new job gets started.
	 * 
	 */
	public static void resetDoneCount () {
		Worker.done = 0;
	}
	
	public Worker (Task task, int totalTasks) {
		super();
		this.task = task;
		Worker.totalTasks = totalTasks;
	}
	
	public String call() {
		if(Parameters.verbose) {
			System.out.println(task.getTaskInfo());
		}
		
		// strand auto detection
		if(this.task.type == Constants.TYPE_STRAND_DETECTION_TASK) {
			StrandDetection.runDetection(task);
		} 
		// else: general scan / target modes
		else {
			if(Parameters.mode.equalsIgnoreCase(Constants.MODE_SCAN)) {
				ScanModeRun.runScanMode(task);
			} else if(Parameters.mode.equalsIgnoreCase(Constants.MODE_TARGET)) {
				TargetModeRun.runTargetMode(task);
			}
		}
		
		
		// call once
		this.task.peakMemory = CheckMemory.checkUsedMemoryMB();
		countTasks();
		
		return task.getTaskInfo();
	}
	
	private static synchronized void countTasks () {
		done++;
		System.out.println(done+"/"+totalTasks);
	}
	
	
}
