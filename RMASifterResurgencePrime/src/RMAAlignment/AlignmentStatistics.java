package RMAAlignment;


import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
/**
 * This class is used to compute some Statistics from a List of Alignments. Alignments should come from CompositionMap and must have their duplicates marked
 * @author huebler
 *
 */
public class AlignmentStatistics {
	private ArrayList<Alignment> currentList;
	private ArrayList<Double> generalStatistics;
	private HashMap<Integer,Integer> coverageHistogram;
	private int length;
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
	public int getLength(){
		return this.length;
	}
	private ArrayList<Alignment> removeDuplicates(ArrayList<Alignment> input){
		HashSet<Integer> length = new HashSet<Integer>();
		if(input != null && input.size() > 2){
			ArrayList<Alignment> positionsToKeep = new ArrayList<Alignment>();
			for(Alignment al : input){
				length.add(al.getReferenceLength());
				if(!al.isDuplicate())
					positionsToKeep.add(al);
			}
			int temp = 0;
			for(int l:length)
				temp+=l;
			this.length = temp;
			return positionsToKeep;
			}
		else return input;
	}
	
	// process best list of start positions
	public void calculateStatistics(){
		
		ArrayList<Alignment> input = removeDuplicates(currentList);
		if(input != null && input.size()>0){
			ArrayList<Double> results = new ArrayList<Double>();
			ReferenceMap refMap = new ReferenceMap(input);
			refMap.process();
			HashMap<Integer,Integer> coverageHistogram = refMap.getCoverageHistogram();
			int temp = 0;
			for(int key : coverageHistogram.keySet())
				temp += coverageHistogram.get(key);
			int zeros = refMap.getLength() - temp ;	
			if(zeros<0)
				zeros = 0;
			coverageHistogram.put(0, zeros);	
			int unique=coverageHistogram.get(1);
			this.coverageHistogram = coverageHistogram;
			results.add(unique/(refMap.getPossible()));
			int nonDuplicates = input.size();
		
			results.add((double) input.size()-refMap.getStackedReads());
			results.add((double) nonDuplicates);
			results.add((double) currentList.size());
			results.add((double) refMap.getLength());
			results.add(refMap.getAverageCoverage());
			this.generalStatistics = results;
		}else{
			HashMap<Integer,Integer> coverageHistogram = new HashMap<Integer,Integer>();
			for(int l = 0; l<=11; l++)
				coverageHistogram.put(l, 0);	
			this.coverageHistogram = coverageHistogram;
			ArrayList<Double> results =	new ArrayList<Double>();
			for(int i = 0; i<5;  i++)
				results.add(0.0);
			this.generalStatistics = results;
		}
		
	}
}
