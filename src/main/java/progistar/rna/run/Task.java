package progistar.rna.run;

import java.util.ArrayList;

import progistar.rna.data.FileInfo;
import progistar.rna.data.Record;

public class Task {

	public ArrayList<Record> records;
	public FileInfo fileInfo;
	
	public Task(ArrayList<Record> records, FileInfo fileInfo) {
		this.records = records;
		this.fileInfo = fileInfo;
	}
}
