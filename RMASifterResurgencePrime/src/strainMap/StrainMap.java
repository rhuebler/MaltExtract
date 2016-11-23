package strainMap;
import java.util.ArrayList;
import java.util.HashMap;

import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;

public class StrainMap {
	/**
	 * Containter of StrainMisMatchContainer and other information for one strain 
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
	public String getEditDistanceHistogram(){
		HashMap<Integer,Integer> histo = new HashMap<Integer,Integer> ();
		for(int i = 0;i < 7;i++){
			histo.put(i, 0);
		}
		ArrayList<Integer> list = container.getEditDistances();
		if(list != null){
			for(int d : list){
				if(d<=5){
					int value = histo.get(d);
					value++;
					histo.put(d, value);
				}else{
					int value = histo.get(6);
					value++;
					histo.put(6, value);
				}	
			}
		}
		
		return name.replace(" ", "_")+"\t"+ histo.get(0)+"\t"+histo.get(1)+"\t"+histo.get(2)+"\t"+histo.get(3)+"\t"+histo.get(4)+"\t"+histo.get(5)+"\t"+histo.get(6);
	}
	public String getPercentIdentityHistogram(){
		HashMap<Integer,Integer> histo = new HashMap<Integer,Integer> ();
		for(int i = 0;i < 5; i++){
			histo.put(i, 0);
		}
		ArrayList<Double> list = container.getPercentIdentity();
		if(list != null){
			for(double d: list){
				if(77.5 <= d && d< 82.5){
					int value = histo.get(0);
					value++;
					histo.put(0, value);
				}else if(82.5 <= d && d< 87.5){
					int value = histo.get(1);
					value++;
					histo.put(1, value);
				}else if(87.5 <= d && d< 92.5){
					int value = histo.get(2);
					value++;
					histo.put(2, value);
				}else if(92.5 <= d && d< 97.5){
					int value = histo.get(3);
					value++;
					histo.put(3, value);
				}else if(97.5 <= d){
					int value = histo.get(4);
					value++;
					histo.put(4, value);
				}
			}
		}
		return name.replace(" ", "_")+"\t"+histo.get(0)+"\t"+histo.get(1)+"\t"+histo.get(2)+"\t"+histo.get(3)+"\t"+histo.get(4);
	}
	public String getReadLengthDistribution(){
		ArrayList<Integer> lengths = container.getLengths();
		if(lengths.size() != 0){
			DescriptiveStatistics stats = new DescriptiveStatistics();
			for(int i : lengths)
				stats.addValue(i);
			String line = name.replace(" ", "_");
			line += "\t"+stats.getMean();
			line += "\t"+stats.getGeometricMean();
			line += "\t"+stats.getPercentile(50);
			line += "\t"+stats.getStandardDeviation();
			return line;
		}else{
			return name.replace(" ", "_")+"\t0\t0\t0\t0";
		}
	}
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
			if(container.getProcessed()!=0){
				container.processMisMatches();
				HashMap<Integer,Double> damage = container.getDamage(); 
				HashMap<Integer,Double> noise =container.getNoise();
		
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
