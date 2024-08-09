package progistar.scan.function;

import java.io.File;
import java.util.Collection;

import org.ahocorasick.trie.Emit;
import org.ahocorasick.trie.Trie;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import progistar.scan.data.Constants;
import progistar.scan.data.LocationInformation;
import progistar.scan.run.Scan;
import progistar.scan.run.Task;

public class ScanModeRun extends Mode{
	
	public static void runScanMode (Task task) {
		if(task.type == Constants.TYPE_SCAN_MODE_TASK) {
			if(Scan.verbose) {
				System.out.println(task.chrName+":"+task.start+"-"+task.end);
			}
			scanReads(task);
		}
	}
	
	private static void scanReads (Task task) {
		long startTime = System.currentTimeMillis();
		// to prevent racing
		File file = new File(Scan.bamFile.getAbsolutePath());
		try (SamReader samReader = SamReaderFactory.makeDefault().open(file)) {
			SAMRecordIterator iterator = null;
			if(task.readType == Constants.MAPPED_READS) {
				iterator = samReader.queryOverlapping(task.chrName, task.start, task.end);
			} else {
				iterator = samReader.queryUnmapped();
			}
			find(iterator, Task.allTrie, task);
		} catch(Exception e) {
			e.printStackTrace();
			System.exit(1);
		}
		
		long endTime = System.currentTimeMillis();
		if(Scan.verbose) {
			System.out.println("Task"+task.taskIdx+" "+(endTime-startTime)/1000+" sec");
		}
	}
}











