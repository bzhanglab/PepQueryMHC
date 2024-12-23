package progistar.scan.parser;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.Hashtable;

import progistar.scan.data.Constants;
import progistar.scan.data.Exon;
import progistar.scan.data.Gene;
import progistar.scan.data.GeneArray;
import progistar.scan.data.Transcript;
import progistar.scan.function.IndexConvertor;

public class ParseGTF {

	// prevent to generate constructor
	private ParseGTF () {}
	
	private static final String[] SELECTED_FEATURES = {"exon", "cds"};
	
	
	private enum FieldIndex {
		CHR(0), SOURCE(1), FEATURE(2), START(3), END(4), STRAND(6), ATTR(8);

		private int value;

		FieldIndex(int value) {
			this.value = value;
		}
	}
	
	private static String getGtfAttr(String[] attr, String tag){
		
		for(String _s : attr){
			if(_s.contains(tag)){
				return _s.replaceAll("[\"\\s]|"+tag, "");
			}
		}
		
		return null;
	}
	
	public static GeneArray[] parseGTF (File gtfFile) {
		System.out.println("Build an interval tree: "+gtfFile.getAbsolutePath());
		long startTime = System.currentTimeMillis();
		Hashtable<String, Gene> geneTable = new Hashtable<String, Gene>();
		Hashtable<String, Transcript> transcriptTable = new Hashtable<String, Transcript>();
		try {
			BufferedReader BR = new BufferedReader(new FileReader(gtfFile));
			
			String line = null;
			
			int chrIndex = FieldIndex.CHR.value;
			int startIndex = FieldIndex.START.value;
			int endIndex = FieldIndex.END.value;
			int featureIndex = FieldIndex.FEATURE.value;
			int strandIndex = FieldIndex.STRAND.value;
			int attrIndex = FieldIndex.ATTR.value;
			
			// add transcript information and corresponding exon structures
			while((line = BR.readLine()) != null) {
				if(line.startsWith("#")) continue; // skip meta
				
				String[] fields = line.split("\t");
				
				String feature = fields[featureIndex];
				int start = Integer.parseInt(fields[startIndex]);
				int end = Integer.parseInt(fields[endIndex]);
				String[] attr = fields[attrIndex].split(";");
				
				if(feature.equalsIgnoreCase("transcript")) {
					String geneID = getGtfAttr(attr, "gene_id");
					String geneName = getGtfAttr(attr, "gene_name");
					String geneType = getGtfAttr(attr, "gene_type");
					String transcriptID = getGtfAttr(attr, "transcript_id");
					
					boolean strand = fields[strandIndex].equalsIgnoreCase("-") ? false : true;
					
					String chr = fields[chrIndex];
					// enroll chr index and chr string
					IndexConvertor.putChrIndexer(chr);
					int chrIdx = IndexConvertor.chrToIndex(chr);
					
					Gene gene = geneTable.get(geneID);
					if(gene == null) {
						gene = new Gene();
						gene.chrIdx = chrIdx;
						gene.start = start;
						gene.end = end;
						gene.id = geneID;
						gene.name = geneName;
						gene.type = geneType;
						geneTable.put(geneID, gene);
					} else {
						gene.start = Integer.min(gene.start, start);
						gene.end = Integer.max(gene.end, end);
					}
					
					Transcript transcript = transcriptTable.get(transcriptID);
					if(transcript == null) {
						transcript = new Transcript();
						transcript.id = transcriptID;
						transcript.start = start;
						transcript.end = end;
						transcript.strand = strand;
						gene.transcripts.put(transcriptID, transcript);
						transcriptTable.put(transcriptID, transcript);
					} else {
						System.out.println(transcriptID +" is duplicated and ignored.");
					}
				} 

				else {
					boolean isSelectedFeature = false;
					for(String sFeature : SELECTED_FEATURES) {
						if(feature.equalsIgnoreCase(sFeature)) isSelectedFeature = true;
					}
					
					// skip if it is not a selected feature.
					if(!isSelectedFeature) continue;
					
					// selected features consist of structural blocks which are building block of a transcript.
					byte bFeature = feature.equalsIgnoreCase("CDS") ? Constants.CDS : Constants.NCDS;
					
					
					String transcriptID = getGtfAttr(attr, "transcript_id");
					Transcript transcript = transcriptTable.get(transcriptID);
					
					// ASSERT!
					Exon aBlock = new Exon();
					aBlock.start = start;
					aBlock.end = end;
					aBlock.feature = bFeature;
					
					transcript.exons.add(aBlock);
				}
			}
			
			BR.close();
			
			// classify exon into CDS, NCDS and UTR.
			// plus, define intron and intergenic regions.
			
		}catch(IOException ioe) {
			System.out.println("...\tFail to load GTF");
		}
		
		

		System.out.println("The number of chromosomes: "+IndexConvertor.size());
		System.out.println("The number of genes: "+geneTable.size());
		System.out.println("The number of transcripts: "+transcriptTable.size());
		GeneArray[] geneArrays = new GeneArray[IndexConvertor.size()];
		geneTable.forEach((id, g)->{
			int chrIdx = g.chrIdx;
			if(geneArrays[chrIdx] == null) {
				geneArrays[chrIdx] = new GeneArray();
				geneArrays[chrIdx].chrIdx = chrIdx;
			}
			geneArrays[chrIdx].genes.add(g);
		});
		
		for(GeneArray geneArray : geneArrays) {
			geneArray.refine();
			for(int i=0; i< geneArray.genes.size(); i++) {
				Gene gene = geneArray.genes.get(i);
				gene.transcripts.forEach((id, t)->{
					t.refine();
				});
			}
		}
		
		
		long endTime = System.currentTimeMillis();
		System.out.println("\tElapsed time: "+((endTime-startTime)/1000) + " sec");
		
		return geneArrays;
		
	}
}
