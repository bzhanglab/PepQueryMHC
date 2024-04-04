package progistar.rna.data;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Hashtable;

import progistar.rna.run.Task;

public class Table {

	private Hashtable<String, Hashtable<String, Integer>> frCnts = new Hashtable<String, Hashtable<String, Integer>>();
	private Hashtable<String, Hashtable<String, Integer>> rcCnts = new Hashtable<String, Hashtable<String, Integer>>();
	private ArrayList<String> fileNames = new ArrayList<String>();
	private ArrayList<String> sampleNames = new ArrayList<String>();
	private Hashtable<String, Double> totalReadCounts = new Hashtable<String, Double>();
	
	public void addStat (Task task) {
		String fileName = task.fileInfo.file.getName();
		// enroll
		// 1) total reads
		// 2) file name / 3) sample name
		totalReadCounts.put(fileName, task.fileInfo.totalReads);
		fileNames.add(fileName);
		sampleNames.add(task.fileInfo.sampleName);
		
		for(Record record : task.records) {
			Hashtable<String, Integer> frMap = frCnts.get(record.frSequence);
			if(frMap == null) {
				frMap = new Hashtable<String, Integer>();
				frCnts.put(record.frSequence, frMap);
			}
			frMap.put(fileName, record.frCnt);			
			
			Hashtable<String, Integer> rcMap = rcCnts.get(record.rcSequence);
			if(rcMap == null) {
				rcMap = new Hashtable<String, Integer>();
				rcCnts.put(record.rcSequence, rcMap);
			}
			rcMap.put(fileName, record.rcCnt);
		}
	}
	
	public void write (File file, ArrayList<Record> records) throws IOException {
		BufferedWriter BW = new BufferedWriter(new FileWriter(file));

		BW.append(Record.header);
		for(String sampleNameName : sampleNames) {
			BW.append("\t").append(sampleNameName);
		}
		
		for(String sampleNameName : sampleNames) {
			BW.append("\tRPHM_").append(sampleNameName);
		}
		
		/*
		 * @Deprecated
		 * 
		for(String fileName : fileNames) {
			BW.append("\tfr+rc_").append(fileName);
		}*/
		BW.newLine();
		
		for(Record record : records) {
			BW.append(record.record);
			Hashtable<String, Integer> frMap = frCnts.get(record.frSequence);
			Hashtable<String, Integer> rcMap = rcCnts.get(record.rcSequence);
			
			for(String fileName : fileNames) {
				int frCnt = frMap.get(fileName);
				int rcCnt = rcMap.get(fileName);
				
				// for unmapped reads
				if(record.chr.equalsIgnoreCase("*")) {
					BW.append("\t"+(frCnt+rcCnt));
				}
				// for mapped reads
				else {
					BW.append("\t"+(frCnt));
				}
			}
			
			for(String fileName : fileNames) {
				int frCnt = frMap.get(fileName);
				int rcCnt = rcMap.get(fileName);
				double totalReads = totalReadCounts.get(fileName);
				/// RPHM
				// for unmapped reads
				if(record.chr.equalsIgnoreCase("*")) {
					BW.append("\t"+Math.pow(10,8) * ((frCnt+rcCnt)/totalReads));
				}
				// for mapped reads
				else {
					BW.append("\t"+Math.pow(10,8) * ((frCnt)/totalReads));
				}
			}
			
			BW.newLine();
		}
		
		BW.close();
	}
}
