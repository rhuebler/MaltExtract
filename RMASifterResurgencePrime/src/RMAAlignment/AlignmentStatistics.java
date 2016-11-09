package RMAAlignment;


import java.util.ArrayList;
import java.util.HashMap;
/**
 * This class is used to compute some Statistics from a List of Alignments. Alignments should come from CompositionMap and must have their duplicates marked
 * @author huebler
 *
 */
public class AlignmentStatistics {
	private ArrayList<Alignment> currentList;
	private ArrayList<Double> generalStatistics;
	private HashMap<Integer,Integer> coverageHistogram;
	public AlignmentStatistics(ArrayList<Alignment> list){
		this.currentList = list;
	}
	//getters
	public ArrayList<Double> getGenaralStatistics(){
		return this.generalStatistics;
	}
	public HashMap<Integer,Integer> getConverageHistogram(){
		return this.coverageHistogram;
	}
	private ArrayList<Alignment> removeDuplicates(ArrayList<Alignment> input){
		if(input != null && input.size() > 2){
			ArrayList<Alignment> positionsToKeep = new ArrayList<Alignment>();
			for(Alignment al : input){
				if(!al.isDuplicate())
					positionsToKeep.add(al);
			}
			return positionsToKeep;
			}
		else return input;
		
	}
	
	// process best list of start positions
	public void calculateStatistics(){
		HashMap<Integer,Integer> coverageHistogram = new HashMap<Integer,Integer>();
		for(int l = 0; l<=11; l++)
			coverageHistogram.put(l, 0);
		ArrayList<Alignment> input = removeDuplicates(currentList);
		if(input != null&&input.size()>2){
			ArrayList<Double> results = new ArrayList<Double>();
			int i = 0; // better solution implemented
			double unique = 0;
			double possible = 0;
			HashMap<Integer,Integer> coverageContainer = new HashMap<Integer,Integer>();
			while(i<input.size()){
				//calculate unique positions per read plus average distance between current and next read
				Alignment current = input.get(i);
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
			
			while(i<input.size()){
			//calculate unique positions per read plus average distance between current and next read
				Alignment current = input.get(i);
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
		int zeros = currentList.get(0).getReferenceLength()-coverageContainer.size();	
		if(zeros<0)
			zeros = 0;
		coverageHistogram.put(0, zeros);	
		unique=coverageHistogram.get(1);
		this.coverageHistogram = coverageHistogram;
		results.add(unique/(possible));
		int nonDuplicates = input.size();
		input.removeAll(stacked);
		results.add((double) input.size());
		results.add((double) nonDuplicates);
		results.add((double) currentList.size());
		results.add((double) currentList.get(0).getReferenceLength());
		this.generalStatistics = results;
		}else{
		this.coverageHistogram = coverageHistogram;
		ArrayList<Double> results =	new ArrayList<Double>();
		for(int i = 0; i<5;  i++)
			results.add(0.0);
		this.generalStatistics = results;
		}
		
	}
}
