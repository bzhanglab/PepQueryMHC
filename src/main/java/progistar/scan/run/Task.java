package progistar.scan.run;

import java.util.ArrayList;

import progistar.scan.data.FileInfo;
import progistar.scan.data.Record;

public class Task {

	public ArrayList<Record> records;
	public FileInfo fileInfo;
	
	public Task(ArrayList<Record> records, FileInfo fileInfo) {
		this.records = records;
		this.fileInfo = fileInfo;
	}
}
