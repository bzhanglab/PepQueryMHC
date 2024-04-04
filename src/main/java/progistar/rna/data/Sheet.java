package progistar.rna.data;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;

public class Sheet {
	
	public ArrayList<FileInfo> fileInfos = new ArrayList<FileInfo>();

	public Sheet(File file) throws IOException {
		System.out.println("Read file information: "+file.getName());
		BufferedReader BR = new BufferedReader(new FileReader(file));
		String line = null;
		
		// skip header
		BR.readLine();
		while((line = BR.readLine()) != null) {
			String[] fields = line.split("\t");
			String sampleName = fields[0];
			String filePath = fields[1];
			String totalReads = fields[2];
			FileInfo fileInfo = new FileInfo();
			
			fileInfo.sampleName = sampleName;
			fileInfo.file = new File(filePath);
			
			if(totalReads.length() != 0) {
				fileInfo.totalReads = Double.parseDouble(totalReads);
			}
			
			fileInfos.add(fileInfo);
		}
		
		BR.close();
	}
}
