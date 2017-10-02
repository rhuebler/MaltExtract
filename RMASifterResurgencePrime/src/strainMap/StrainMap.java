package strainMap;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.HashMap;
import RMAAlignment.AlignmentStatistics;

public class StrainMap {
	/**
	 * Containter of StrainMisMatchContainer and other information for one strain 
	 * 
	 */
	//Initialize attributes
	private String name;
	private StrainMisMatchContainer container;
	private int numMatches=0;
	private  AlignmentStatistics statistics;
	
	//constructor
	public StrainMap(String s,StrainMisMatchContainer container, int n){
		this.name = s;
		this.container = container;
		this.numMatches = n;
		
	}
	//setters
	//get edit edit distance and process edidistances
	public void setStatistics(){
		 AlignmentStatistics stats = container.getStatistics();
		 this.statistics = stats;
		
	}
	public String getCoverageHistogram(){
		HashMap<Integer,Integer>histogram = statistics.getConverageHistogram();
		String coveragHistograms =  name +"\t" + name;
		for(int k : histogram.keySet())
			coveragHistograms += "\t" + histogram.get(k);
		return coveragHistograms;
	}
	public String getCoveragePositions(){
		String coveragePostitions =  name +"\t" + name;
		for(String cov :statistics.getCoveragePositions()){
			coveragePostitions += "\t"+cov;
		}
		return coveragePostitions;
	}
	public String getReadDistribution(){
		DecimalFormat df = new DecimalFormat("#.###");
		String readDist = name +"\t" + name;
		for(double d :statistics.getGenaralStatistics())
			readDist += "\t" + df.format(d);
		return readDist;
	}
	public String getEditDistanceHistogram(){
		HashMap<Integer,Integer> histo = new HashMap<Integer,Integer> ();
		for(int i = 0;i < 11;i++){
			histo.put(i, 0);
		}
		ArrayList<Integer> list = container.getEditDistances();
		if(list != null){
			for(int d : list){
				if(d<=10){
					int value = histo.get(d);
					value++;
					histo.put(d, value);
				}else{
					int value = histo.get(11);
					value++;
					histo.put(11, value);
				}	
			}
		}
		
		return name.replace(" ", "_")+"\t"+ histo.get(0)+"\t"+histo.get(1)+"\t"+histo.get(2)+"\t"+histo.get(3)+"\t"+histo.get(4)+
				"\t"+histo.get(5)+"\t"+histo.get(6)+"\t"+histo.get(7)+"\t"+histo.get(8)+"\t"+histo.get(9)+"\t"+histo.get(10)+"\t"+histo.get(11);
	}
	// get and process percent identity
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
	private int round(double i, int v){
	    return (int) (Math.round(i/v) * v);
	}
	// produce and get read distances
	public String getReadLengthDistribution(){
		ArrayList<Integer> lengths = container.getLengths();
		if(lengths.size() != 0){
			HashMap<Integer,Integer> intervals = new HashMap<Integer,Integer>();
			for(int i = 25;i<=200;i+=5)
				intervals.put(i, 0);
				for(int i : lengths){// round to the closest number dividable by 5 to get the intervals
					int value = intervals.get(round(i, 5));
					value++;
					intervals.replace(round(i, 5), value);
				}
				String rlDist = name;
			for(int key:intervals.keySet()){
					rlDist+="\t"+intervals.get(key);
				}
			return rlDist;
		}else{
			return name.replace(" ", "_")+"\t0\t0\t0\t0";
		}
	}
	//get strain mismatch container
	public StrainMisMatchContainer getStrainMisMatchContainer(){
		return this.container;
	}
	public String getName(){
		return this.name;
	}
	public int getNumberOfMatches(){
		return this.numMatches;
	}
	public String getDamageLine(){ // process Map Into Damage Output Line 
			if(container.getProcessed()!=0){
				container.processMisMatches();
				numMatches = container.getProcessed();
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
	//setters
	public void setStrainMisMatchContainer(StrainMisMatchContainer container){
		this.container = container;
	}
	public void setNumberOfMatches(int n){
		this.numMatches = n;
	}
	
	
	
}
