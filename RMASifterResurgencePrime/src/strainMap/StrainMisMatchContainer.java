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
	String name;


	private int processed = 0;
	// constructor
	public  StrainMisMatchContainer() {
	}
	public void setName(String name){
		this.name =  name;
	}
//getters
	public int getProcessed(){
		return this.processed;
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
public String getDamageLine(){ // process Map Into Damage Output Line 
	if(getProcessed()!=0){
		processMisMatches();
		int numMatches = getProcessed();
		HashMap<Integer,Double> damage = getDamage(); 
		HashMap<Integer,Double> noise = getNoise();

		String part1 = name.replace(" ", "_");
		String part2 = "";
		for(int i = 0;i < 20; i++){
			if(damage.containsKey(i)){
				part1 += "\t" + damage.get(i);
			}else{	
				part1 += "\t" + 0;	
			}
			if(noise.containsKey(i)){
				part2 += "\t" + noise.get(i);
			}else{
				part2 += "\t" + 0;
			}	
		}
		if(numMatches!=0)
			part1 += part2 + "\t" + numMatches;
		else
			part1 += part2 + "\t" + 0;
		return part1;
}else{
	String part1 = name.replace(" ", "_");
	for(int i = 0;i<=40;i++){
		part1+="\t"+0;
	}
	return part1;
}
}
}
