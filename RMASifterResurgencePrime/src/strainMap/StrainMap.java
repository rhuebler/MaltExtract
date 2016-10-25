package strainMap;
import java.util.HashMap;

public class StrainMap {
	/**
	 * containter of StrainMisMatchContainer for one strain
	 */
	private String name;
	private StrainMisMatchContainer container;
	private int numMatches;
	
	//constructor
	public StrainMap(String s,StrainMisMatchContainer container, int n){
		this.name = s;
		this.container = container;
		this.numMatches = n;
		
	}
	//setters
	public void setStrainMisMatchContainer(StrainMisMatchContainer container){
		this.container = container;
	}
	public void setNumberOfMatches(int n){
		this.numMatches = n;
	}
	
	// getters
	public StrainMisMatchContainer getStrainMisMatchContainer(){
		return this.container;
	}
	public String getName(){
		return this.name;
	}
	public int getNumberOfMatches(){
		return this.numMatches;
	}
	public String getLine(){ // process Map Into Damage Output Line
		container.processMisMatches();
		HashMap<Integer,Double> damage = container.getDamage(); 
		HashMap<Integer,Double> noise =container.getNoise();
		
		String part1 = name;
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
	}
}
