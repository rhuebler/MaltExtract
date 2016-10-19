package strainMap;
// Store all Combinations for the first 10 and last 10 Positions of each Blast hit
import java.util.HashMap;

import RMAAlignment.Alignment;

public class StrainMisMatchContainer {
	private HashMap<Integer,HashMap<String,Integer>> container = new HashMap<Integer,HashMap<String,Integer>>();
public void process(Alignment al){
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
}
