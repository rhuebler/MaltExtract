package RMAAlignment;


import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
/**
 * This class is used to mark and remove stacked reads from alignments. Alignments should come from CompositionMap and must have their duplicates marked. 
 * Return list with removed stacked reads and coverage distribution. Is an version of alignment statistics that does not retrieve statistics but instead 
 * just removed overlapping reads
 * @author huebler
 *
 */
public class GetStackedReads {
	// intilaize with values
	private ArrayList<Alignment> currentList;
	private ArrayList<Double> generalStatistics;
	private HashMap<Integer,Integer> coverageHistogram;
	private int length;
	private ArrayList<Alignment> nonStacked = new ArrayList<Alignment>();
	private boolean turnedOn = true;
	private boolean turnOffDeDupping = false;
	private boolean turnOffDestacking = false;

	public GetStackedReads(boolean turnOffDeDupping, boolean turnOffDestacking){
		this.turnOffDeDupping = turnOffDeDupping;
		this.turnOffDestacking = turnOffDestacking;
	}
	//getters
	public boolean wasTurnedOn(){
		return this.turnedOn;
	}
	public  GetStackedReads(ArrayList<Alignment> list){
		this.currentList = list;
	}
	public ArrayList<Alignment> getNonStacked(){
		return this.nonStacked;
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
	//remove duplicate alignments
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
			input = removeDuplicates(currentList);
		else
			 input = removeDuplicates(currentList);
		if(input != null && input.size()>0){
			ReferenceMap refMap = new ReferenceMap(input,turnOffDestacking);
			refMap.process();
			this.nonStacked = refMap.getNonStacked();
			this.turnedOn = refMap.wasTurnedOn();
		}
		
	}
}
