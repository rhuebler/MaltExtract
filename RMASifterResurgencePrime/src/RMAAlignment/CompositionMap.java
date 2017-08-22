package RMAAlignment;
import java.util.ArrayList;
import java.util.HashMap;
/**
 * Save for a Megan Node all Blast hits stored as Alignment Objects these are kept as hashmap
 *  of the ID of the matching species and a List of the Hits 
 *  Mainly serves as a Storage to identify the reference to calculate Read distribution and to mark duplicates 
 */
// set Composition Map
public class CompositionMap {
	private boolean turnOffDestacking= false;
	private boolean turnOffDeDupping = false;
//	public CompositionMap(HashMap<Integer,ArrayList<Alignment>> map){
//		setCompositionMap(map);
//	}
	public CompositionMap(HashMap<Integer, HashMap<String, ArrayList<Alignment>>> taxonMap,boolean turnOffDestacking, boolean turnOffDeDupping){
		setCompositionMap(taxonMap);
		this.turnOffDestacking = turnOffDestacking;
		this.turnOffDeDupping = turnOffDeDupping;
	}
	// initialize values
private HashMap<Integer, HashMap<String, ArrayList<Alignment>>> compositionMap;// hashMap of ReferenceID to List of start positions
private HashMap<String,ArrayList<Alignment>> resultsMap = new HashMap<String,ArrayList<Alignment>>();// hashMap of ReferenceID to List of start positions
private int maxID;
private String maxReference;
private ArrayList<Double> generalStatistics;
private HashMap<Integer,Integer> coverageHistogram;
private ArrayList<String> coveragePositions;
private int refLength = 0;
private boolean turnedOn = true;
//getter
public ArrayList<String> getCoveragePositions(){
	return this.coveragePositions;
}
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

public HashMap<Integer, HashMap<String, ArrayList<Alignment>>> getCompositionMap(){
	return this.compositionMap;}
public ArrayList<Double> getGenaralStatistics(){
	return this.generalStatistics;
}
public HashMap<Integer,Integer> getConverageHistogram(){
	return this.coverageHistogram;
}
//setters
private void setCompositionMap(HashMap<Integer, HashMap<String, ArrayList<Alignment>>> taxonMap){// technically not required
	this.compositionMap = taxonMap;}

private void setMaxID(int i){
	this.maxID = i;
}
private void setMaxReference(String reference){
	this.maxReference = reference;
}
//utility
private HashMap<String,ArrayList<Alignment>> markDuplicates(HashMap<String,ArrayList<Alignment>> map){ // that seems to be correct may be the ordering works better and therefore more stacking reads can be found
	AlignmentComparator comp = new AlignmentComparator(); // should theoretically sort my alignments according to start positions
	for(String key : map.keySet()){
		ArrayList<Alignment> inList = map.get(key);
		inList.sort(comp);
		int i =0;
		HashMap<String, Integer> timesViewed =  new HashMap<String, Integer>();
		while(i < inList.size() - 1){// array size is 47 last 46
			
			Alignment current = inList.get(i);
			if(timesViewed.containsKey(current.getQuery())){
				int k=timesViewed.get(current.getQuery());
				timesViewed.replace(current.getQuery(),k );
			}else{
				timesViewed.put(current.getQuery(), 1);
			}
			
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
			if( cStart == nStart && cEnd == nEnd && current.getReferenceLength() == next.getReferenceLength()){
				if(timesViewed.containsKey(next.getQuery())||current.getQuery().equals(next.getQuery())){
				inList.get(i+1).setDuplicate(true);
				}
			}
				i++;
		}
		map.replace(key, inList);
	}
	return map;
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

	AlignmentStatistics stats = new AlignmentStatistics(compositionMap.get(getMaxID()).get(maxReference),turnOffDestacking, turnOffDeDupping);
	stats.calculateStatistics();
	this.coverageHistogram = stats.getConverageHistogram();
	this.generalStatistics = stats.getGenaralStatistics();
	this.refLength = stats.getLength();
	this.coveragePositions = stats.getCoveragePositions();
}

// process composition and find taxon with maximum number of start positions
public void getNonStacked(){
	//get nonstacked Reads
	for(int key : compositionMap.keySet()){
		HashMap<String,ArrayList<Alignment>> rMap = compositionMap.get(key);
		for(String reference : rMap.keySet()){
			if(turnOffDestacking){
				for(Alignment al : rMap.get(reference)){
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
			GetStackedReads reads = new GetStackedReads(rMap.get(reference));
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
}


// process compostion map
public void process(){
	HashMap<Integer,HashMap<String,ArrayList<Alignment>>> cMap = compositionMap;
	int max=0; 
	int maxKey=0; 
	String maxReference="";
	for(int id: cMap.keySet()){
		HashMap<String,ArrayList<Alignment>> map = cMap.get(id);
		for(String key : map.keySet()){
			if(map.get(key).size()>max){
				maxKey=id;
				max= map.get(key).size();
				maxReference=key;
			}
		}
	}
	setMaxID(maxKey);
	setMaxReference(maxReference);
	markAllDuplicates();
  }

//calculate and get all top references
public HashMap<Double,ArrayList<Integer>> getAllTopReferences(){
	HashMap<Double,ArrayList<Integer>> additional = new HashMap<Double,ArrayList<Integer>>();
	HashMap<Integer,HashMap<String,ArrayList<Alignment>>> map = compositionMap;
	if(map.size()>=1){
		int maxSize = map.get(maxID).get(maxReference).size();
		for(double x = 1.0;x>=0.0;x-=0.1){
			ArrayList<Integer> margin = new ArrayList<Integer>();
			for(int key : map.keySet()){
				int size = 0;
				for(String ref : map.get(key).keySet()){
					size += map.get(key).get(ref).size();
				}		
				if(size > (maxSize*x) && size <= maxSize*(x+0.1)){
				margin.add(key);
				}	
			}
			additional.put(x, margin);
		}	
	}
	return additional;
}

}//class
