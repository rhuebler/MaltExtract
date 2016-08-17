package RMAAlignment;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
/**
 * Save for a Megan Node all Blast hits stored as Alignment Objects these are kept as hashmap
 *  of the ID of the matching species and a List of the Hits 
 */
public class CompositionMap {
	public CompositionMap(HashMap<Integer,ArrayList<Alignment>> map){
		setCompositionMap(map);
	}
private HashMap<Integer,ArrayList<Alignment>> compositionMap;// hashMap of ReferenceID to List of start positions
//HashMap<Integer,Integer> taxonComposition; currently unused 
private int maxID;
//from stack overflow
private double getMean(ArrayList<Integer> data)
{
    double sum = 0.0;
    for(int a : data)
        sum += a;
    return (sum/data.size());
}
//TODO add ALexanders moving window function 
private double getVariance(ArrayList<Integer> data, double mean){
    double temp = 0;
    for(double a :data)
        temp += (mean-a)*(mean-a);
    return temp/data.size();}

private double getStdDev(double variance)
{
    return Math.sqrt(variance);
}

private double getMedian(ArrayList<Integer> data) 
{	 double pos1 = Math.floor((data.size() - 1.0) / 2.0);
     double pos2 = Math.ceil((data.size() - 1.0) / 2.0);
     if (pos1 == pos2 ) {
        return data.get((int)pos1);
     } else {
        return (data.get((int)pos1) + data.get((int)pos2)) / 2.0 ;
     }
}
private ArrayList<OrientedRead> removeDuplicates(ArrayList<OrientedRead> inList){// seems like relatively good and clean way to solve it 
	int i =0;//bit annoying how often i reiterate through the same lists maybe this can be retconned 
	ArrayList<OrientedRead> posToRemove = new ArrayList<OrientedRead> ();
	while(i < inList.size() - 1){// array size is 47 last 46
		if(inList.get(i).getStart() == inList.get(i+1).getStart()
		&& inList.get(i).getEnd() == inList.get(i+1).getEnd()){
			posToRemove.add(inList.get(i));
			posToRemove.add(inList.get(i+1));
		}
			i++;
		
			
			inList.removeAll(posToRemove);
	}
	return inList;
}
//owncode
public int getMaxID(){
	return this.maxID;}

public HashMap<Integer,ArrayList<Alignment>> getCompositionMap(){
	return this.compositionMap;}

private void setCompositionMap(HashMap<Integer,ArrayList<Alignment>> map){// technically not rrequired
	this.compositionMap = map;}

private void setMaxID(int i){
	this.maxID = i;
}
// process best list of start positions
public List<Double> getStatistics(){
	ArrayList<Alignment> input = getCompositionMap().get(getMaxID());
	ArrayList<OrientedRead> list = new ArrayList<OrientedRead>();
	int refLength = input.get(0).getReferenceLength();
	ArrayList<Integer> distance = new ArrayList<Integer>();
	List<Double> results= new ArrayList<Double>();
	for(Alignment al : input)
		list.add(new OrientedRead(al.getStart(),al.getEnd()));
	
	ReadComparator comp = new ReadComparator(); // should theoretically sort my alignments according to start positions
	list.sort(comp);// great works
	list = removeDuplicates(list);
	int i = 0;//TODO maybe there is a better solution also is there an logical error in here 
	double unique = 0;
	while(i<list.size()){
		OrientedRead current = list.get(i);
		int length = current.getEnd() - current.getStart() + 1;
		int soli = 0;
		if(i==0){
			OrientedRead next = list.get(i+1);
			int a = (current.getEnd()-next.getStart());
			if(a > 0)
				soli = length - a;
			else if(a == 0)
				soli = length - 1;
			else 
				soli = length;
			distance.add(next.getStart()-current.getStart());
			}
		else if(i<list.size()-1){
			OrientedRead next = list.get(i+1);
			OrientedRead previous = list.get(i-1);
			int a = current.getEnd() - next.getStart(); 
			int b =	previous.getEnd() - current.getStart();
			if(a<0 && b<0)
				soli = length;
			else if(a<0 && b>0)
				soli = length - a;
			else if(a>0 && b<0)
				soli = length - b;
			if(a==0)soli -= 1;// test if either distance is zero and one position less can be considered unique 
			if(b==0)soli -= 1;
			distance.add(next.getStart()-current.getStart());
		}
		else{
			OrientedRead previous = list.get(i-1);
			int a = previous.getEnd() - current.getStart();
			if(a>0)
				soli = length-a;
				else if(a==0)
					soli = length-1;
				else if(a<0)
					soli= length;
		}
			
		if(soli <= 0 || soli >= length)
			soli = length;
		unique+=soli;
		i++;
	}
	
	double mean = getMean(distance);
	double median = getMedian(distance);
	double variance = getVariance(distance,mean);
	double std = getStdDev(variance);
	results.add(mean);
	results.add(median);
	results.add(variance);
	results.add(std);
	results.add((unique/refLength)/(100.0*list.size())); //TODO consider taking actual average match length
	results.add((double)list.size());
	results.add((double) input.size());
	results.add((double)input.get(0).getReferenceLength());
	return results;
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
	//setTaxonCompositionMap(results);
  }

//private void setTaxonCompositionMap(HashMap<Integer, Integer> results) {
//	this.taxonComposition = results;
//	
//}	
//public HashMap<Integer,Integer> getTaxonComposition(){
//	return this.taxonComposition;
//}

}//class
// class that can represent an alingment and switch it around internally depending on orientation 
class OrientedRead{
	private int start;
	private int end;
	public OrientedRead(int start, int end){
		this.end = end;
		this.start = start;
		checkOrientation();
	}
	public int getStart(){
		return this.start;
	}
	public int getEnd(){
		return this.end;
		}
	private void checkOrientation() // should I use in other classes as well 
	{ 
	 int start = getStart();
	 int end = getEnd();
	 if(start>end){
		this.start = end;
		this.end = start;		
	   }
		
	}
}

class ReadComparator implements Comparator<OrientedRead> // comparator for OrientedRead class
{
  @Override public int compare( OrientedRead re1, OrientedRead re2 ) // allows to sort them by reordered starting positions 
  {
    return re1.getStart()-re2.getStart();
  }
}