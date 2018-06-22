package RMAAlignment;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
/**
 * This class is used to get the coverage Histogram and to mark and remove stacked reads based on input parameters
 * Essenitally this is used to "map" reads and get the coverage and find out which reads overlap
 * @author huebler
 *
 */
public class ReferenceMap {
	//initialize attributes
	private ArrayList<Alignment> input;
	private double coverageDeviation=0.0;
	private HashMap<Integer,Integer> coverageHistogram;
	private HashMap<Integer,Double> coverageMap = new HashMap<Integer,Double>();
	private int stackedReads = 0;
	private double unique = 0;
	private double possible = 0;
	private int length = 0;
	private ArrayList<Alignment> nonStacked = new ArrayList<Alignment>();
	private boolean turnedOn = true;
	private double averageCoverage = 0;
	private boolean turnOffDestacking = false;
	// constructor
	public ReferenceMap(ArrayList<Alignment> input, boolean turnOffDestacking){
		this.input = input;
		this.turnOffDestacking = turnOffDestacking;
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
	//differentiate different plasmids and chromosomes by length and calculate total length  while processing the alignments
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
		double averageCoverage = 0;//Get information on coverage
		for(i=0;i<5;i++){
			coverageMap.put(i,0.0);
		}	
		// get percentage higher than x
		for(int k: coverageContainer.keySet()){
			double coverage = coverageContainer.get(k);
			averageCoverage+=coverage;
			if(coverage>=1){
				Double count = coverageMap.get(0);
				count++;
				coverageMap.replace(0, count);
			}if(coverage>=2){
				Double count = coverageMap.get(1);
				count++;
				coverageMap.replace(1, count);
			}if(coverage>=3){
				Double count = coverageMap.get(2);
				count++;
				coverageMap.replace(2, count);
			}if(coverage>=4){
				Double count = coverageMap.get(3);
				count++;
				coverageMap.replace(3, count);
			}if(coverage>=5){
				Double count = coverageMap.get(4);
				count++;
				coverageMap.replace(4, count);
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
		if(averageCoverage>=5 && !turnOffDestacking)
			turnedOn = false;
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
			l++;
		}
		stackedReads+=stacked.size();
		nonStacked.removeAll(stacked);
	}

}
