package strainMap;
// Store all Combinations for the first 10 and last 10 Positions of each Blast hit
import java.util.HashMap;

import RMAAlignment.Alignment;

public class StrainMisMatchContainer {
	private HashMap<Integer,HashMap<String,Integer>> container = new HashMap<Integer,HashMap<String,Integer>>();
public void processAlignment(Alignment al){
	if(al.getMlength() >= 20){
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
					System.out.println(s);
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
					System.out.println(s);
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
public void processMisMatches(){
	HashMap<Integer,Double> damage = new HashMap<Integer,Double>();
	HashMap<Integer,Double> noise = new HashMap<Integer,Double>();
	for(int i = 0; i<20;i++){
		HashMap<String, Integer>positionContainer = container.get(i);
		if(i<10){
			int damDivident = 0;
			int damDivisor = 0;
			if(positionContainer.containsKey("CT")){
				damDivisor=positionContainer.get("CT");
				if(positionContainer.containsKey("CC")){
					damDivident = (positionContainer.get("CT")+positionContainer.get("CC"));
					double d = (damDivisor/damDivident);
					damage.put(i, d);
				}else{
					double d = 1.0;
					damage.put(i, d);
				}
			}	
			int divisor = 0;
			double divident = 0.0;
			for(String key : positionContainer.keySet()){
				if(!key.equals("CT")){
					divisor+=positionContainer.get(key);
				}
			}
			divident= divisor;
			if(positionContainer.containsKey("GG")){
				divident-=positionContainer.get("GG");
			}
			if(positionContainer.containsKey("AA")){
				divident-=positionContainer.get("AA");
			}
			if(positionContainer.containsKey("TT")){
				divident-=positionContainer.get("TT");
			}
			divident -= damDivident;
			noise.put(i, divident/divisor);
		}else{
			int damDivident = 0;
			int damDivisor = 0;
			if(positionContainer.containsKey("GA")){
				damDivident=positionContainer.get("GA");
				if(positionContainer.containsKey("GG")){
					damDivisor = (positionContainer.get("GA")+positionContainer.get("GG"));
					double d = (damDivident/damDivisor);
					damage.put(i, d);
				}else{
					double d = 1.0;
					damage.put(i, d);
				}
			}
			double divident = 0.0;
			int divisor = 0;
			for(String key : positionContainer.keySet()){
				if(!key.equals("GA"))
					divisor+=positionContainer.get(key);
			}
			divident= divisor;
			if(positionContainer.containsKey("CC")){
				divident-=positionContainer.get("CC");
			}
			if(positionContainer.containsKey("AA")){
				divident-=positionContainer.get("AA");
			}
			if(positionContainer.containsKey("TT")){
				divident-=positionContainer.get("TT");
			}
			divident -= damDivident;
			noise.put(i, divident/divisor);
			
		}
	}
}
}
