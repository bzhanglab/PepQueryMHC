package progistar.scan.function;

import java.io.File;
import java.io.IOException;

import htsjdk.samtools.util.BlockCompressedInputStream;

public class ReadBGZF {

	public static void main(String[] args) throws IOException {
		BlockCompressedInputStream reader = new BlockCompressedInputStream(new File("/Users/seunghyukchoi/eclipse-workspace/BAMScanner/test.bgzf"));
		String line = null;
		
		int count = 10;
		while((line = reader.readLine()) != null) {
			System.out.println(line);
			
			if(--count == 0) {
				break;
			}
		}
		
	}
}
