package RMAAlignment;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;

public class ReferenceMap {
	private ArrayList<Alignment> input;
	private HashMap<Integer,Integer> coverageHistogram;
	private int stackedReads = 0;
	private double unique = 0;
	private double possible = 0;
	private int length = 0;
	private ArrayList<Alignment> nonStacked = new ArrayList<Alignment>();
	private boolean turnedOn = true;
	private double averageCoverage = 0;
	public ReferenceMap(ArrayList<Alignment> input){
		this.input = input;
	}
	//getters
	public double getAverageCoverage(){
		return averageCoverage;
	}
	public boolean wasTurnedOn(){
		return this.turnedOn;
	}
	public ArrayList<Alignment> getNonStacked(){
		return this.nonStacked;
	}
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
	//differentiate different plasmids and chromosomes by length and calculate total length 
	public void process(){
		nonStacked.addAll(input);
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
			calculate(totalReferences.get(key),key);
			length+=key;
		}		
		
	}
	
	private void calculate (ArrayList<Alignment> alignments, int length){
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
		HashSet<Alignment> stacked = new HashSet<Alignment>();
		
		//calculate average coverage on reference and turnOff stackedReadRemoval if too high
		double averageCoverage = 0;
		for(int k: coverageContainer.keySet())
			averageCoverage+=coverageContainer.get(k);
		averageCoverage/=length;
		this.averageCoverage = averageCoverage;
		if(averageCoverage>=10)
			turnedOn = false;
		
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
				if(coverage>1 && turnedOn){
					current.setStacked(true);
					stacked.add(current);
				}
				if(coverage <= 10){
					coverageHistogram.replace(coverage, coverageHistogram.get(coverage)+1);	
				}else if(coverage >= 11){
					coverageHistogram.replace(11, coverageHistogram.get(11)+1);
					
				}
			}
			i++;
		}
		stackedReads+=stacked.size();
		nonStacked.removeAll(stacked);
	}

}
