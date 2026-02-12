package progistar.scan.function;

import java.io.File;
import java.io.IOException;

import htsjdk.samtools.BAMIndexer;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

public class BAMUtil {

	public static void index (File file) throws IOException {
		File indexFile = new File(file.getAbsolutePath()+".bai");
		
		if(indexFile.exists()) {
			return;
		}
		
		try (SamReader reader = SamReaderFactory.makeDefault().enable(SamReaderFactory.Option.INCLUDE_SOURCE_IN_RECORDS).open(file)) {
            if (!reader.type().equals(SamReader.Type.BAM_TYPE)) {
                throw new IllegalArgumentException("Input file is not a BAM file");
            }
            System.out.println("Indexing BAM file...");
            BAMIndexer indexer = new BAMIndexer(indexFile, reader.getFileHeader());
            reader.iterator().forEachRemaining(indexer::processAlignment);
            indexer.finish();

            System.out.println("Index created: " + indexFile.getAbsolutePath());
        } catch (Exception e) {
            e.printStackTrace();
        }
	}
	
	public static File sort (File file) throws IOException {
		File outputFile = new File(file.getAbsolutePath().replace(".bam", ".sorted.bam"));
		
		System.out.println("Sort BAM file...");
		// 1. 입력 파일 열기
        try (SamReader reader = SamReaderFactory.makeDefault().open(file)) {
            SAMFileHeader header = reader.getFileHeader();
            
            // 2. 정렬 순서를 'coordinate'로 설정 (인덱싱 필수 조건)
            header.setSortOrder(SAMFileHeader.SortOrder.coordinate);

            // 3. 최적화된 Writer 설정
            SAMFileWriterFactory factory = new SAMFileWriterFactory()
                .setMaxRecordsInRam(500000);          // 메모리 내 보유 레코드 수 조절

            try (SAMFileWriter writer = factory.makeBAMWriter(header, false, outputFile)) {
                for (SAMRecord record : reader) {
                    writer.addAlignment(record);
                }
            }
        } catch (Exception e) {
            e.printStackTrace();
        }
        
        System.out.println("Sorted BAM file was created: "+outputFile.getName());
        System.out.println("Delete unsorted BAM...");
        file.delete();
        
        return outputFile;
	}
}
