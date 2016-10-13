package strainMap;
import java.util.HashMap;

public class StrainMap {
	/**
	 * containter of Strain Map and Mismatches
	 */
	private String name;
	private HashMap<Integer,Integer> misMap;
	private HashMap<Integer,Integer> substitutionMap;
	private int numMatches;
	
	//constructor
	public StrainMap(String s,HashMap<Integer,Integer> m,HashMap<Integer,Integer> sub, int n){
		this.name = s;
		this.misMap = m;
		this.numMatches = n;
		this.substitutionMap = sub;
		
	}
	//setters
	public void setMisMap(HashMap<Integer,Integer> m){
		this.misMap = m;
	}
	public void setSubstitutionMap(HashMap<Integer,Integer> m){
		this.substitutionMap = m;
	}
	public void setNumberOfMatches(int n){
		this.numMatches = n;
	}
	
	// getters
	public HashMap<Integer,Integer> getMisMap(){
		return this.misMap;
	}
	public HashMap<Integer,Integer> getSubstitutionMap(){
		return this.substitutionMap;
	}
	public String getName(){
		return this.name;
	}
	public int getNumberOfMatches(){
		return this.numMatches;
	}
	public String getLine(){ // process Map Into Damage Output Line
		String part1 = name;
		String part2 = "";
		for(int i = 0;i < 20; i++){
			if(misMap.containsKey(i)){
				part1 += "\t" + misMap.get(i);
			}else{	
				part1 += "\t" + 0;	
			}
			if(substitutionMap.containsKey(i)){
				part2 += "\t" + substitutionMap.get(i);
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
