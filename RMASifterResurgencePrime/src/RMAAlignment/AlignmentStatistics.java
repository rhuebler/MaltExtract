package RMAAlignment;


import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
/**
 * This class is used to compute some Statistics from a List of Alignments. Alignments should come from CompositionMap and must have their duplicates marked
 * Usually duplicates will be removed unless the off switch is turned on. And stacking reads is removed also retrieve information whether filter was turned off or not
 * @author huebler
 * @params ArrayList<Alignment> list, boolean turnOffDestacking, boolean turnOffDeDupping
 */
public class AlignmentStatistics {
	// initialize attributes
	private ArrayList<Alignment> currentList;
	private ArrayList<Double> generalStatistics;
	private HashMap<Integer,Integer> coverageHistogram;
	private int length;
	private boolean turnOffDestacking = false;
	private boolean turnOffDeDupping = false;
	private ArrayList<String> coverage = new ArrayList<String>();
	private ArrayList<Alignment> destackedList;
	//constructor and set values
	public AlignmentStatistics(ArrayList<Alignment> list, boolean turnOffDestacking, boolean turnOffDeDupping){
		this.currentList = list;
		this.turnOffDestacking = turnOffDestacking;
		this.turnOffDeDupping = turnOffDeDupping;
	}
	//getters
	public ArrayList<Alignment> getDestackedList(){
		return destackedList;
	}
	public ArrayList<String> getCoveragePositions(){
		return this.coverage;
	}
	public ArrayList<Double> getGenaralStatistics(){
		return this.generalStatistics;
	}
	public HashMap<Integer,Integer> getConverageHistogram(){
		return this.coverageHistogram;
	}
	public int getLength(){
		return this.length;
	}
	// remove reads that are flagged as duplicates
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
		
		ArrayList<Alignment> input = new ArrayList<Alignment>();
		if(!turnOffDeDupping)
			input =removeDuplicates(currentList);
		else
			input = currentList;
		if(input != null && input.size()>0){
			ArrayList<Double> results = new ArrayList<Double>();
			// initialize reference map
			ReferenceMap refMap = new ReferenceMap(input,turnOffDestacking);
			refMap.process();
			HashMap<Integer,Integer> coverageHistogram = refMap.getCoverageHistogram();
			int temp = 0;
			for(int key : coverageHistogram.keySet())
				temp += coverageHistogram.get(key);
			int zeros = refMap.getLength() - temp ;	
			if(zeros<0)
				zeros = 0;
			
			DecimalFormatSymbols otherSymbols = new DecimalFormatSymbols();
			otherSymbols.setDecimalSeparator('.');
			otherSymbols.setGroupingSeparator(','); 
			DecimalFormat df = new DecimalFormat("#.###",otherSymbols);
			coverageHistogram.put(0, zeros);
			int unique=coverageHistogram.get(1);
			
			if(unique  == 0){
				results.add(0.0);
			}else{
				results.add((Double.parseDouble(df.format(unique/(refMap.getPossible())))));
			}
			this.coverageHistogram = coverageHistogram;
			destackedList = refMap.getNonStacked();
			results.add((double) input.size()-refMap.getStackedReads());
			results.add((double)  input.size());
			results.add((double) currentList.size());
			results.add((double) refMap.getLength());

			
			coverage.add(df.format(refMap.getAverageCoverage()));
			coverage.add(df.format(refMap.getCoverageDeviation()));
			for(int key:refMap.getCoveragePositions().keySet()){
				coverage.add(df.format(refMap.getCoveragePositions().get(key)/refMap.getLength()));
			}
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
