package RMAAlignment;
import java.util.ArrayList;
import java.util.HashMap;
/**
 * Save for a Megan Node all Blast hits stored as Alignment Objects these are kept as hashmap
 *  of the ID of the matching species and a List of the Hits 
 *  Mainly serves as a Storage to identify the reference to calculate Read distribution and to mark duplicates 
 */
public class CompositionMap {
	private boolean turnOffDestacking= false;
	public CompositionMap(HashMap<Integer,ArrayList<Alignment>> map){
		setCompositionMap(map);
	}
	public CompositionMap(HashMap<Integer,ArrayList<Alignment>> map,boolean turnOffDestacking){
		setCompositionMap(map);
		this.turnOffDestacking = turnOffDestacking;
	}
private HashMap<Integer,ArrayList<Alignment>> compositionMap;// hashMap of ReferenceID to List of start positions
private HashMap<String,ArrayList<Alignment>> resultsMap = new HashMap<String,ArrayList<Alignment>>();// hashMap of ReferenceID to List of start positions
//HashMap<Integer,Integer> taxonComposition; currently unused 
private int maxID;
private ArrayList<Double> generalStatistics;
private HashMap<Integer,Integer> coverageHistogram;
private int refLength = 0;
private boolean turnedOn = true;
//getter
public boolean wasTurnedOn(){
	return this.turnedOn;
}
public HashMap<String,ArrayList<Alignment>> getResultsMap(){
	return this.resultsMap;
}
public int getReferenceLength(){
	return this.refLength;
}
public int getMaxID(){
	return this.maxID;}

public HashMap<Integer,ArrayList<Alignment>> getCompositionMap(){
	return this.compositionMap;}
public ArrayList<Double> getGenaralStatistics(){
	return this.generalStatistics;
}
public HashMap<Integer,Integer> getConverageHistogram(){
	return this.coverageHistogram;
}
//setters
private void setCompositionMap(HashMap<Integer,ArrayList<Alignment>> map){// technically not required
	this.compositionMap = map;}

private void setMaxID(int i){
	this.maxID = i;
}
//utility
private ArrayList<Alignment> markDuplicates(ArrayList<Alignment> inList){ // that seems to be correct may be the ordering works better and therefore more stacking reads can be found
	AlignmentComparator comp = new AlignmentComparator(); // should theoretically sort my alignments according to start positions
	inList.sort(comp);
	int i =0;
	while(i < inList.size() - 1){// array size is 47 last 46
		Alignment current = inList.get(i);
		Alignment next = inList.get(i+1);
		int cStart = 0;
		int cEnd = 0;
		int nStart = 0;
		int nEnd = 0;
		if(!current.isReversed() == !next.isReversed()){
			cStart = current.getStart();
			cEnd = current.getEnd();
			nStart = next.getStart();
			nEnd = next.getEnd();
		}else if(current.isReversed() == !next.isReversed()){
			cStart = current.getEnd();
			cEnd = current.getStart();
			nStart = next.getStart();
			nEnd = next.getEnd();
		}else if(!current.isReversed() == next.isReversed()){
			cStart = current.getStart();
			cEnd = current.getEnd();
			nStart = next.getEnd();
			nEnd = next.getStart();
		}else{
			cStart = current.getEnd();
			cEnd = current.getStart();
			nStart = next.getEnd();
			nEnd = next.getStart();
		}
		if( cStart == nStart && cEnd == nEnd &&
		 current.getReferenceLength() == next.getReferenceLength()){
			//inList.get(i).setDuplicate(true);TODO only set first Read as duplicate 
			inList.get(i+1).setDuplicate(true);
		}
			i++;
	}
	return inList;
}
public void markReferenceDuplicates(){
	this.compositionMap.put(getMaxID(), markDuplicates(this.compositionMap.get(getMaxID())));
}
public void markAllDuplicates(){
	for(int key : this.compositionMap.keySet()){
	this.compositionMap.put(key, markDuplicates(this.compositionMap.get(key)));
	}
}
public void calculateStatistics(){
	AlignmentStatistics stats = new AlignmentStatistics(this.compositionMap.get(getMaxID()));
	stats.calculateStatistics();
	this.coverageHistogram = stats.getConverageHistogram();
	this.generalStatistics = stats.getGenaralStatistics();
	this.refLength = stats.getLength();
}

// process composition and find taxon with maximum number of start positions
public void getNonStacked(){
	//get nonstacked Reads
	for(int key : compositionMap.keySet()){
		if(turnOffDestacking){
			for(Alignment al : compositionMap.get(key)){
				String name = al.getReadName();
				if(resultsMap.containsKey(name)){
					ArrayList<Alignment> list = resultsMap.get(name);
					list.add(al);
					resultsMap.replace(name, list);
				}else{
					ArrayList<Alignment> list = new ArrayList<Alignment>();
					list.add(al);
					resultsMap.put(name, list);
					}
			}
		}else{
			GetStackedReads reads = new GetStackedReads(compositionMap.get(key));
			reads.calculateStatistics();
			if(!reads.wasTurnedOn())
				turnedOn = false;
			for(Alignment al : reads.getNonStacked()){
				String name = al.getReadName();
				if(resultsMap.containsKey(name)){
					ArrayList<Alignment> list = resultsMap.get(name);
					list.add(al);
					resultsMap.replace(name, list);
				}else{
					ArrayList<Alignment> list = new ArrayList<Alignment>();
					list.add(al);
					resultsMap.put(name, list);
					}
				}
			}
		}
	}

public void process(){
	HashMap<Integer,ArrayList<Alignment>> map = compositionMap;
	HashMap<Integer,Integer> results = new HashMap<Integer,Integer>();
	int max=0; int maxKey=0;
	for(int key : map.keySet()){
		results.put(key, map.get(key).size());
		if(map.get(key).size()>max){
			maxKey=key;
			max= map.get(key).size();
			}
	}
	setMaxID(maxKey);
	markAllDuplicates();
	//setTaxonCompositionMap(results);
  }
public ArrayList<Integer> getAllTopReferences(){
	ArrayList<Integer> additional = new ArrayList<Integer>();
	HashMap<Integer,ArrayList<Alignment>> map = compositionMap;
	if(map.size()>1){
		int maxSize = map.get(maxID).size();
		for(int key : map.keySet()){
			if(key != getMaxID() && map.get(key).size()>= (maxSize-maxSize*0.25)){
				additional.add(key);
			}
		}
	}
	return additional;
}

}//class
