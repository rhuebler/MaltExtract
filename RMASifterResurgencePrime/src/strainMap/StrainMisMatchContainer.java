package strainMap;
import java.util.ArrayList;
// Store all Combinations for the first 10 and last 10 Positions of each Blast hit
import java.util.HashMap;

import RMAAlignment.Alignment;
import RMAAlignment.AlignmentComparator;
import RMAAlignment.AlignmentStatistics;
import behaviour.Filter;
/**
 * Functions that collects all misMatches and processes them 
 * also serves as a container for alignmets
 * @author huebler
 *
 */

public class StrainMisMatchContainer{
	// initialize attributes
	private HashMap<Integer,HashMap<String,Integer>> container = new HashMap<Integer,HashMap<String,Integer>>();
	private HashMap<Integer,Double> damage; 
	private HashMap<Integer,Double> noise;
	
	private ArrayList<Integer> distances;
	private ArrayList<Double> pIdents;
	private ArrayList<Integer> lengths;
	private ArrayList<Alignment> alignments;

	private int processed = 0;
	private Filter filter;
	// constructor
	public  StrainMisMatchContainer(Filter filter) {
		this.filter = filter;
		if(filter == Filter.CRAWL){
			distances = new ArrayList<Integer>();
			pIdents = new ArrayList<Double>();
			lengths = new ArrayList<Integer>();
			alignments = new ArrayList<Alignment>();
		}
	}
	private ArrayList<Alignment> markDuplicates(ArrayList<Alignment>list){ // that seems to be correct may be the ordering works better and therefore more stacking reads can be found
		AlignmentComparator comp = new AlignmentComparator(); // should theoretically sort my alignments according to start positions
			ArrayList<Alignment> inList = list;
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
		return inList;
	}
//getters
	public int getProcessed(){
		return this.processed;
	}
	public ArrayList<Integer> getEditDistances(){
		return this.distances;
	}
	public ArrayList<Integer> getLengths(){
		return this.lengths;
	}
	public ArrayList<Double> getPercentIdentity(){
		return this.pIdents;
	}
	public AlignmentStatistics getStatistics(){
		AlignmentStatistics stats = new AlignmentStatistics(markDuplicates(alignments),false,false);
		stats.calculateStatistics();
		return stats;
	}
public HashMap<Integer,Double> getDamage(){
	return this.damage;
}
public HashMap<Integer,Double> getNoise(){
	return this.noise;
}	
//retrieve information from alignment
public void processAlignment(Alignment al){
	if(al.getMlength() >= 20){
		if(filter == Filter.CRAWL){
		distances.add(al.getEditDistance());
		pIdents.add(al.getPIdent());
		lengths.add(al.getMlength());
		alignments.add(al);
		}
		
		processed+=1;
		String q = al.getQuery();
		String r = al.getReference();
		for(int i = 0; i< 20; i++){
			if(i<10){
			if(container.containsKey(i)){
					HashMap<String,Integer> misMatches = container.get(i);
					String s =	r.charAt(i)+""+q.charAt(i);
					if(misMatches.containsKey(s)){
					misMatches.replace(s, misMatches.get(s)+1);
					}else{
						misMatches.put(s, 1);
					}
				}else{
					HashMap<String,Integer> misMatches = new HashMap<String,Integer>();
					String s =	r.charAt(i)+""+q.charAt(i);
					misMatches.put(s, 1);
					container.put(i, misMatches);
				}
			}else{
				if(container.containsKey(i)){
					HashMap<String,Integer> misMatches = container.get(i);
					String s =	r.charAt(al.getMlength()+i-20)+""+q.charAt(al.getMlength()+i-20);
					if(misMatches.containsKey(s)){
					misMatches.replace(s, misMatches.get(s)+1);
					}else{
						misMatches.put(s, 1);
					}
					//System.out.println(s);
				}else{
					HashMap<String,Integer> misMatches = new HashMap<String,Integer>();
					String s =	r.charAt(i)+""+q.charAt(i);
					misMatches.put(s, 1);
					container.put(i, misMatches);
				}
			}
		}
	}
}
// get mismatches
public void processMisMatches(){
	HashMap<Integer,Double> damage = new HashMap<Integer,Double>();
	HashMap<Integer,Double> noise = new HashMap<Integer,Double>();
	for(int i = 0; i<20;i++){
		HashMap<String, Integer>positionContainer = container.get(i);
		if(i<10){
			double damDivident = 0.0;
			double damDivisor = 0.0;
			if(positionContainer.containsKey("CT")){
				damDivisor=positionContainer.get("CT");
				if(positionContainer.containsKey("CC")){
				
					damDivident = (damDivisor+positionContainer.get("CC"));
					
					double d = damDivisor/damDivident;
					damage.put(i, d);
				}else{
					double d = 1.0;
					damage.put(i, d);
				}
			}	
			double divisor = 0.0;
			double divident = 0.0;
			for(String key : positionContainer.keySet()){
				if(!key.equals("CT") |! key.contains("-")){
					divisor+=positionContainer.get(key);
				}
			}
			divident += divisor;
			if(positionContainer.containsKey("GG")){
				divident-=positionContainer.get("GG");
			}
			if(positionContainer.containsKey("AA")){
				divident-=positionContainer.get("AA");
			}
			if(positionContainer.containsKey("TT")){
				divident-=positionContainer.get("TT");
			}
			if(positionContainer.containsKey("CC")){
				divident-=positionContainer.get("CC");
			}
			divident /= 11;
			noise.put(i, divident/divisor);
		}else{
			double damDivident = 0.0;
			double damDivisor = 0.0;
			if(positionContainer.containsKey("GA")){
				damDivident=positionContainer.get("GA");
				if(positionContainer.containsKey("GG")){
					damDivisor = (damDivident+positionContainer.get("GG"));
					double d = damDivident/damDivisor;
					damage.put(i, d);
				}else{
					double d = 1.0;
					damage.put(i, d);
				}
			}
			double divident = 0.0;
			double divisor = 0.0;
			for(String key : positionContainer.keySet()){
				if(!key.equals("GA") || !key.contains("-"))
					divisor+=positionContainer.get(key);
			}
			divident += divisor;
			if(positionContainer.containsKey("CC")){
				divident-=positionContainer.get("CC");
			}
			if(positionContainer.containsKey("AA")){
				divident-=positionContainer.get("AA");
			}
			if(positionContainer.containsKey("TT")){
				divident-=positionContainer.get("TT");
			}
			if(positionContainer.containsKey("GG")){
				divident-=positionContainer.get("GG");
			}
			divident /= 11;
			noise.put(i, divident/divisor);
			
		}
	}
	this.damage = damage;
	this.noise = noise;
}
}
