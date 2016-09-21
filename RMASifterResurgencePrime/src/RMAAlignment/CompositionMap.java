package RMAAlignment;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
/**
 * Save for a Megan Node all Blast hits stored as Alignment Objects these are kept as hashmap
 *  of the ID of the matching species and a List of the Hits 
 *  Mainly serves as a Storage to identify the reference to calculate Read distribution and to mark duplicates 
 */
public class CompositionMap {
	public CompositionMap(HashMap<Integer,ArrayList<Alignment>> map){
		setCompositionMap(map);
	}
private HashMap<Integer,ArrayList<Alignment>> compositionMap;// hashMap of ReferenceID to List of start positions
//HashMap<Integer,Integer> taxonComposition; currently unused 
private int maxID;
//from stack overflow


//owncode
public int getMaxID(){
	return this.maxID;}

public HashMap<Integer,ArrayList<Alignment>> getCompositionMap(){
	return this.compositionMap;}

private void setCompositionMap(HashMap<Integer,ArrayList<Alignment>> map){// technically not required
	this.compositionMap = map;}

private void setMaxID(int i){
	this.maxID = i;
}
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
		if( cStart == nStart && cEnd == nEnd){
			inList.get(i).setDuplicate(true);
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
public List<Double> getStatistics(){
	AlignmentStatistics stats = new AlignmentStatistics(this.compositionMap.get(getMaxID()));
	return(stats.getStatistics());
}
// process composition and find taxon with maximum number of start positions
public void process(){
	HashMap<Integer,ArrayList<Alignment>> map = getCompositionMap();
	HashMap<Integer,Integer> results = new HashMap<Integer,Integer>();
	int max=0; int maxKey=0;
	for(int key : map.keySet()){
		results.put(key, map.get(key).size());
		if(map.get(key).size()>=max){
			maxKey=key;
			max= map.get(key).size();
			}
	}
	setMaxID(maxKey);
	markAllDuplicates();
	//setTaxonCompositionMap(results);
  }

}//class
