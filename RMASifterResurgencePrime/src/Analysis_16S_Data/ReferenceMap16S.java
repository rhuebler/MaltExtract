package Analysis_16S_Data;

import java.util.ArrayList;
import java.util.HashMap;

import RMAAlignment.Alignment;

public class ReferenceMap16S {

	//initialize attributes
	private ArrayList<Alignment> input;
	private double coverageDeviation=0.0;
	private HashMap<Integer,Integer> coverageHistogram;
	private HashMap<Integer,Double> coverageMap = new HashMap<Integer,Double>();
	private int stackedReads = 0;
	private int length = 0;
	private double averageCoverage = 0;
	// constructor
	public ReferenceMap16S(ArrayList<Alignment> input){
		this.input = input;
	}
	//getters
	public double getCoverageDeviation(){
		return this.coverageDeviation;
	}
	public HashMap<Integer,Double>  getCoveragePositions(){
		return coverageMap;
	}
	public double getAverageCoverage(){
		return averageCoverage;
	}
	public int getLength(){
		return this.length;
	}
	public int getStackedReads(){
		return this.stackedReads;
	}
	public HashMap<Integer,Integer> getCoverageHistogram(){
		return this.coverageHistogram;
	}
	//differentiate different plasmids and chromosomes by length and calculate total length  while processing the alignments
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
		//calculate average coverage on reference and turnOff stackedReadRemoval if too high
		double averageCoverage = 0;//Get information on coverage
		for(i=0;i<5;i++){
			coverageMap.put(i,0.0);
		}	
		// get percentage higher than x
		for(int k: coverageContainer.keySet()){
			double coverage = coverageContainer.get(k);
			averageCoverage+=coverage;
			if(coverage>=1){
				coverageMap.replace(0,coverageMap.get(0)+1);
			}if(coverage>=2){
				coverageMap.replace(1,coverageMap.get(1)+1);
			}if(coverage>=3){
				coverageMap.replace(2, coverageMap.get(2)+1);
			}if(coverage>=4){
				coverageMap.replace(3, coverageMap.get(3)+1);
			}if(coverage>=5){
				coverageMap.replace(4, coverageMap.get(4)+1);
			}
		}
		averageCoverage/=length;
		// calculate standard deviation
		double temp = 0.0;
		this.averageCoverage = averageCoverage;
		for(int k: coverageContainer.keySet()){
			double coverage = coverageContainer.get(k);
			temp+=(coverage-averageCoverage)*(coverage-averageCoverage);
		}	
		temp+=(averageCoverage*(length-coverageContainer.size()));
		temp/=length;
		this.coverageDeviation = temp;
		int l =0;
		while(l<alignments.size()){
		//calculate unique positions per read plus average distance between current and next read
			Alignment current = alignments.get(l);
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
					int coverage = coverageContainer.get(k);
					if(coverage <= 10){
						coverageHistogram.replace(coverage, coverageHistogram.get(coverage)+1);	
					}else if(coverage >= 11){
						coverageHistogram.replace(11, coverageHistogram.get(11)+1);
					}
			}
			l++;
		}
	}

}
