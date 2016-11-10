package RMAAlignment;

import java.util.ArrayList;
import java.util.HashMap;

public class ReferenceMap {
	private ArrayList<Alignment> input;
	private HashMap<Integer,Integer> coverageHistogram;
	private int stackedReads = 0;
	private double unique = 0;
	private double possible = 0;
	private int length = 0;
	public ReferenceMap(ArrayList<Alignment> input){
		this.input = input;
	}
	//getters
	public int getLength(){
		return this.length;
	}
	public double getUnique(){
		return this.unique;
	}
	public double getPossible(){
		return this.possible;
	}
	public int getStackedReads(){
		return this.stackedReads;
	}
	public HashMap<Integer,Integer> getCoverageHistogram(){
		return this.coverageHistogram;
	}
	public void process(){
		HashMap<Integer,Integer> coverageHistogram = new HashMap<Integer,Integer>();
		for(int l = 0; l<=11; l++)
			coverageHistogram.put(l, 0);
		this.coverageHistogram = coverageHistogram;
		HashMap<Integer,ArrayList<Alignment>> totalReferences = new HashMap<Integer,ArrayList<Alignment>>();
		for (Alignment al : input){
			if(totalReferences.containsKey(al.getReferenceLength())){
				ArrayList<Alignment> alignments = totalReferences.get(al.getReferenceLength());
				alignments.add(al);
				totalReferences.replace(al.getReferenceLength(), alignments);
			}else{
				ArrayList<Alignment> alignments = new ArrayList<Alignment>();
				alignments.add(al);
				totalReferences.put(al.getReferenceLength(), alignments);
			}
		}
		for(int key : totalReferences.keySet()){
			calculate(totalReferences.get(key));
			length+=key;
		}		
		
	}
	
	private void calculate (ArrayList<Alignment> alignments){
		int i = 0; // better solution implemented
		HashMap<Integer,Integer> coverageContainer = new HashMap<Integer,Integer>();
		while(i<alignments.size()){
			//calculate unique positions per read plus average distance between current and next read
			Alignment current = alignments.get(i);
			int cStart = 0;
			int cEnd = 0;
			if(current.isReversed()){
				cStart = current.getEnd();
				cEnd = current.getStart();		
			}else{
				cStart = current.getStart();
				cEnd = current.getEnd();		
			}
			for(int k = cStart; k<= cEnd; k++){
				if(coverageContainer.containsKey(k)){
					coverageContainer.replace(k, coverageContainer.get(k)+1);
				}else{
					coverageContainer.put(k, 1);
				}
			}
			i++;
		}	
		i = 0;
		ArrayList<Alignment> stacked = new ArrayList<Alignment>();
		
		while(i<alignments.size()){
		//calculate unique positions per read plus average distance between current and next read
			Alignment current = alignments.get(i);
			int cStart = 0;
			int cEnd = 0;
			if(current.isReversed()){
				cStart = current.getEnd();
				cEnd = current.getStart();		
			}else{
				cStart = current.getStart();
				cEnd = current.getEnd();		
			}
			possible += (cEnd - cStart)+1;
			for(int k = cStart; k<= cEnd; k++){
				int coverage = coverageContainer.get(k);
				if(coverage == 1){
					coverageHistogram.replace(1, coverageHistogram.get(1)+1);
				}else if(coverage == 2){
					coverageHistogram.replace(2, coverageHistogram.get(2)+1);
					stacked.add(current);
				}else if(coverage == 3){
					coverageHistogram.replace(3, coverageHistogram.get(3)+1);
					stacked.add(current);
				}else if(coverage == 4){
					coverageHistogram.replace(4, coverageHistogram.get(4)+1);
					stacked.add(current);
				}else if(coverage == 5){
					coverageHistogram.replace(5, coverageHistogram.get(5)+1);
					stacked.add(current);
				}else if(coverage == 6){
					coverageHistogram.replace(6, coverageHistogram.get(6)+1);
					stacked.add(current);
				}else if(coverage == 7){
					coverageHistogram.replace(7, coverageHistogram.get(7)+1);
					stacked.add(current);
				}else if(coverage == 8){
					coverageHistogram.replace(8, coverageHistogram.get(8)+1);
					stacked.add(current);
				}else if(coverage == 9){
					coverageHistogram.replace(9, coverageHistogram.get(9)+1);
					stacked.add(current);
				}else if(coverage == 10){
					coverageHistogram.replace(10, coverageHistogram.get(10)+1);
					stacked.add(current);
				}else if(coverage >= 11){
					coverageHistogram.replace(11, coverageHistogram.get(11)+1);
					stacked.add(current);
				}
			}
			i++;
		}
		stackedReads+=stacked.size();
	}

}
