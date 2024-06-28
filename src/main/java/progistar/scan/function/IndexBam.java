package progistar.scan.function;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Hashtable;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.util.BlockCompressedOutputStream;

public class IndexBam {

	public static void main(String[] args) throws IOException {
		long startTime = System.currentTimeMillis();
		// to prevent racing
		File[] files = new File("/Volumes/Zhanglab/2024_LUAD").listFiles();
		
		// check file
		ArrayList<File> nFiles = new ArrayList<File>();
		for(File file : files) {
			if(file.getName().startsWith(".")) continue;
			if(file.getName().endsWith("Aligned.sortedByCoord.out.bam")) {
				nFiles.add(file);
			}
		}
		
		BlockCompressedOutputStream bgzfOut = new BlockCompressedOutputStream(new File("test.bgzf"));
		Hashtable<String, long[]> readCount = new Hashtable<String, long[]>();
		int fIdx = 0;
		for(File file : nFiles) {
			try (SamReader samReader = SamReaderFactory.makeDefault().open(file)) {
				SAMRecordIterator iterator = samReader.queryOverlapping("chr1", 1, 500000000);
				
				int count = 0;
				while (iterator.hasNext()) {
					SAMRecord samRecord = iterator.next();
					if(samRecord.isSecondaryAlignment()) continue;
					count++;
					
					String sequence = samRecord.getReadString();
					String referenceName = samRecord.getReferenceName();
					int start = samRecord.getAlignmentStart();
					int end = samRecord.getAlignmentEnd();
					String cigar = samRecord.getCigarString();
					
					String key = sequence+"\t"+referenceName+"\t"+start+"\t"+end+"\t"+cigar;
					
					long[] reads = readCount.get(key);
					if(reads == null) {
						reads = new long[nFiles.size()];
						readCount.put(key, reads);
					}
					reads[fIdx]++;
				}
				
				System.out.println(count + "=>" +readCount.size());
			} catch(Exception e) {
				e.printStackTrace();
				System.exit(1);
			}
			
			fIdx++;
		}
		
		// header
		
		bgzfOut.write( ("Sequence\tReferenceName\tStart\tEnd\tCigar").getBytes());
		for(File file : nFiles) {
			bgzfOut.write(("\t"+file.getName()).getBytes());
		}
		bgzfOut.write(System.lineSeparator().getBytes());
		
		readCount.forEach((record, values)->{
			try {
				bgzfOut.write(record.getBytes());
				for(long value : values) {
					bgzfOut.write( ("\t"+value).getBytes() );
				}
				bgzfOut.write(System.lineSeparator().getBytes());
			}catch(IOException ioe) {
				
			}
		});
		
		
		bgzfOut.close();
		long endTime = System.currentTimeMillis();
		System.out.println((endTime-startTime)/1000+" sec");
	}
	
	
	// below codes are deprecated
	public static String encoding (String sequence) {
		StringBuilder encode = new StringBuilder();
		
		int len = sequence.length();
		char prev = sequence.charAt(0);
		int cnt = 0;
		for(int i=1; i<len; i++) {
			char cur = sequence.charAt(i);
			if(prev != cur) {
				encode.append(prev);
				if(cnt != 0) {
					encode.append(cnt);
				}
				prev = cur;
				cnt = 0;
			} else {
				cnt++;
			}
		}
		// last pang!
		encode.append(prev);
		if(cnt != 0) {
			encode.append(cnt);
		}
		
		return encode.toString();
	}
	
	public static String decoding (String encode) {
		StringBuilder sequence = new StringBuilder();
		
		int len = encode.length();
		char prev = encode.charAt(0);
		int cnt = 0;
		for(int i=1; i<len; i++) {
			
		}
		// last pang!
		
		
		return encode.toString();
	}
	
	public static String printHeader(SAMFileHeader header) {
        // Print the version, sort order, and group order
		StringBuilder headerStr = new StringBuilder();
		headerStr.append("@HD\tVN:" + header.getVersion() + "\tSO:" + header.getSortOrder());
		headerStr.append(System.lineSeparator());
		
        // Print reference sequences
        header.getSequenceDictionary().getSequences().forEach(sequence -> {
        	headerStr.append("@SQ\tSN:" + sequence.getSequenceName() + "\tLN:" + sequence.getSequenceLength());
        	headerStr.append(System.lineSeparator());
        });

        // Print read groups
        header.getReadGroups().forEach(readGroup -> {
        	headerStr.append("@RG\tID:" + readGroup.getId() + "\tSM:" + readGroup.getSample());
        	headerStr.append(System.lineSeparator());
        });

        // Print program records
        header.getProgramRecords().forEach(program -> {
        	headerStr.append("@PG\tID:" + program.getId() + "\tCL:" + program.getCommandLine());
        	headerStr.append(System.lineSeparator());
        });

        // Print other header lines
        header.getComments().forEach(comment -> {
        	headerStr.append("@CO\t" + comment);
        	headerStr.append(System.lineSeparator());
        });
        
        return headerStr.toString();
    }
	
}
