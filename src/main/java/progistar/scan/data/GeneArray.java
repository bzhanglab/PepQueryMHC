package progistar.scan.data;

import java.util.ArrayList;
import java.util.Collections;
import java.util.LinkedList;

public class GeneArray {
	
	public int chrIdx;
	public ArrayList<Gene> genes = new ArrayList<Gene>();
	public int mark = 0;
	
	public void refine() {
		Collections.sort(genes);
		
		// init minmax
		for(Gene gene : genes) {
			gene.min = gene.start;
			gene.max = gene.end;
		}
		
		int lastIdx = genes.size()-1;
		refine(lastIdx/2, 0, lastIdx);
		mark++;
	}
	
	private int[] refine(int mid, int left, int right) {
		int lastIdx = genes.size()-1;
		int[] minmax = new int[2];
		
		Gene gene = genes.get(mid);
		
		if(gene.mark == mark) {
			minmax[0] = gene.min;
			minmax[1] = gene.max;
			return minmax;
		}
		
		gene.mark = mark;
		
		minmax[0] = gene.min;
		minmax[1] = gene.max;
		
		int lRight = mid -1;
		int lMid = (left+lRight) / 2;
		if(lRight >= 0) {
			int[] lMinmax = refine(lMid, left, lRight);
			minmax[0] = Math.min(lMinmax[0], minmax[0]);
			minmax[1] = Math.max(lMinmax[1], minmax[1]);
		}
		
		int rLeft = mid + 1;
		int rMid = (rLeft+right) / 2;
		if(rLeft <= lastIdx) {
			int[] rMinmax = refine(rMid, rLeft, right);
			minmax[0] = Math.min(minmax[0], rMinmax[0]);
			minmax[1] = Math.max(minmax[1], rMinmax[1]);
		}
		
		
		gene.min = minmax[0];
		gene.max = minmax[1];
		
		return minmax;
	}
	
	public ArrayList<Gene> findOverlap (int start, int end) {
		
		// if mark exceeds the limited number
		// reset mark
		if(mark == Integer.MAX_VALUE) {
			for(Gene gene : this.genes) {
				gene.mark = -1;
			}
			mark = 0;
		}
		
		mark ++;
		ArrayList<Gene> matchedGenes = new ArrayList<Gene>();
		
		
		LinkedList<Integer> nextMid = new LinkedList<Integer>();
		LinkedList<Integer> nextLeft = new LinkedList<Integer>();
		LinkedList<Integer> nextRight = new LinkedList<Integer>();
		int lastIdx = genes.size()-1;
		
		nextMid.add(lastIdx/2);
		nextLeft.add(0);
		nextRight.add(lastIdx);
		
		while(!nextMid.isEmpty()) {
			int mid = nextMid.pollFirst();
			int left = nextLeft.pollFirst();
			int right = nextRight.pollFirst();
			
			// check tie
			int gStart = genes.get(mid).start;
			int minMid = mid;
			int maxMid = mid;
			// left ties
			for(int idx=mid-1; idx>left; idx--) {
				Gene tieGene = genes.get(idx);
				if(tieGene.start == gStart) {
					minMid = idx;
				} else {
					break;
				}
			}
			// right ties
			for(int idx=mid+1; idx<right; idx++) {
				Gene tieGene = genes.get(idx);
				if(tieGene.start == gStart) {
					maxMid = idx;
				} else {
					break;
				}
			}
			
			for(mid=minMid; mid<=maxMid; mid++) {
				// if overlap?
				Gene midGene = genes.get(mid);
				// is visited?
				if(midGene.mark == mark) {
					continue;
				}
				midGene.mark = mark;
				
				// possible?
				if( (start >= midGene.min && start <= midGene.max) ||
						(end >= midGene.min && end <= midGene.max) ||
						(start <= midGene.min && end >= midGene.max) ) {
					//System.out.println(mid+", "+midGene.geneName+", "+midGene.start + ":" + midGene.end +", "+midGene.min+" & "+midGene.max);
					
					// is overlapped?
					if( (start >= midGene.start && start <= midGene.end) ||
						(end >= midGene.start && end <= midGene.end) ||
						(start <= midGene.start && end >= midGene.end) ) {
						matchedGenes.add(midGene);
					}
					
					// add left
					int lRight = mid -1;
					int lMid = (left+lRight) / 2;
					
					if(lRight >= 0) {
						nextMid.add(lMid);
						nextLeft.add(left);
						nextRight.add(lRight);
					}
					
					// add right
					int rLeft = mid + 1;
					int rMid = (rLeft+right) / 2;
					
					if(rLeft <= lastIdx) {
						nextMid.add(rMid);
						nextLeft.add(rLeft);
						nextRight.add(right);
					}
				}
				
			}
			
		}
		return matchedGenes;
	}
	
}
