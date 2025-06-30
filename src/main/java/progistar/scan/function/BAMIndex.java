package progistar.scan.function;

import java.io.File;
import java.io.IOException;

import htsjdk.samtools.BAMIndexer;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

public class BAMIndex {

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
}
