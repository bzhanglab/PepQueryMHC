package progistar.scan.function;

import java.io.File;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import progistar.scan.data.Constants;
import progistar.scan.data.Parameters;
import progistar.scan.run.Task;

public class StrandDetection {

	public static void runDetection (Task task) {
		if(task.type == Constants.TYPE_STRAND_DETECTION_TASK) {
			if(Parameters.verbose) {
				System.out.println(task.chrName+":"+task.start+"-"+task.end);
			}
			detect(task);
		}
	}
	
	private static void detect (Task task) {
		long startTime = System.currentTimeMillis();
		// to prevent racing
		File file = new File(Parameters.bamFile.getAbsolutePath());
		try (SamReader samReader = SamReaderFactory.makeDefault().open(file)) {
			SAMRecordIterator iterator = samReader.queryOverlapping(task.chrName, 1, Integer.MAX_VALUE);
			int size = task.end;
			
			while((size--) > 0 && iterator.hasNext()) {
				SAMRecord samRecord = iterator.next();
				
				boolean isPass = false;
				if(Parameters.count.equalsIgnoreCase(Constants.COUNT_PRIMARY) && samRecord.isSecondaryOrSupplementary()) {
            		isPass = true;
            	}
				
				Object xsTag = samRecord.getAttribute("XS");
				if(xsTag == null) {
					isPass = true;
				}
				
				if(isPass) {
					continue;
				}
				
				
				int flags = samRecord.getFlags();
				boolean isFirstSegment = (0x40 & flags) == 0x40 ? true : false;
				boolean isForward = (0x10 & flags) == 0x10 ? false : true;
				boolean strand = ((Character) xsTag) == '+' ? true : false;
				
				// first segment
				if(isFirstSegment) {
					if(isForward == strand) {
						task.R1F++;
					} else {
						task.R1R++;
					}
				} 
				// second segment
				else {
					if(isForward == strand) {
						task.R2F++;
					} else {
						task.R2R++;
					}
				}
				
			}
			
			iterator.close();
		} catch(Exception e) {
			e.printStackTrace();
			System.exit(1);
		}
		
		long endTime = System.currentTimeMillis();
		if(Parameters.verbose) {
			System.out.println("Task"+task.taskIdx+" "+(endTime-startTime)/1000+" sec");
		}
	}
}
