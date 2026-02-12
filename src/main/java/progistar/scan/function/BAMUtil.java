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
        try (SamReader reader = SamReaderFactory.makeDefault().open(file)) {
            SAMFileHeader header = reader.getFileHeader();
            
            header.setSortOrder(SAMFileHeader.SortOrder.coordinate);

            SAMFileWriterFactory factory = new SAMFileWriterFactory()
                .setMaxRecordsInRam(500000);

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
