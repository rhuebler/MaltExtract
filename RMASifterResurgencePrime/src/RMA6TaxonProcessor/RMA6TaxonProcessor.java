package RMA6TaxonProcessor;

import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.logging.Level;
import java.util.logging.Logger;

import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;

import NCBI_MapReader.NCBI_MapReader;
import RMAAlignment.Alignment;
import RMAAlignment.CompositionMap;
import behaviour.Filter;
import megan.data.IMatchBlock;
import strainMap.StrainMisMatchContainer;

/**
 * Extract all information from one Taxon and save the information in the specified slots to be retrieved
 * in RMA6Processor is parent class for all filtering Taxonprocessors now mainly serves as parent class for
 * all Filters so that all of them have access to the same basic functions and can be used somewhat interchangeably
 * Initlazed with some standard values.
 * @author huebler
 *
 */
public class RMA6TaxonProcessor {
	/**
	 * @param int ID, NCBI_MapReader reader, boolean verbose, Logger, log, Logger warning, minPIdent
	 * @return int numMatches, String readDistribution, HashMap EditDistance, HashMap Percent Identity
	 */ 
protected String taxName="unasigned_Name";	
protected String readDistribution;
protected NCBI_MapReader mapReader;
protected Integer taxID;
protected double minPIdent;
protected double topPercent;
protected int maxLength;
protected HashMap<Integer,Integer> editHistogram;
protected HashMap<Integer,Integer> pIdentHistogram;
protected boolean verbose;
protected Logger log;
protected Logger warning;
protected ArrayList<String> alignments = new ArrayList<String>();
protected ArrayList<String> alignmentList;
protected ArrayList<String> readList;
protected ArrayList<String> lines = new ArrayList<String>();
protected String damageLine;
protected String coverageLine;
protected int numOfReads = 0;
protected int numMatches = 0;
protected int originalNumberOfReads = 0;
protected int originalNumberOfAlignments = 0;
protected String filterLine = "";
protected ArrayList<Integer> distances = new ArrayList<Integer>();
protected ArrayList<Double> pIdents = new ArrayList<Double>();
protected HashMap<Integer,HashMap<String, ArrayList<Alignment>>> taxonMap = new HashMap<Integer,HashMap<String, ArrayList<Alignment>>>();
protected Filter filter = Filter.NON;
protected StrainMisMatchContainer container = new StrainMisMatchContainer();
protected String readLengthDistribution;
protected String readLengthStatistics;
protected int refLength = 0;
protected ArrayList<Integer> lengths = new ArrayList<Integer>();
protected boolean turnedOn = true;
protected String additionalEntries;
protected String covPositions;
protected DecimalFormat df;
//constructor
public RMA6TaxonProcessor(Integer id, double pID, NCBI_MapReader reader, boolean v, Logger log, Logger warning, double topPercent, int maxLength, Filter f){
	// set input values
	this.mapReader = reader;
	this.minPIdent = pID;
	this.taxID = id;
	this.verbose = v;
	this.log = log;
	this.warning = warning;
	this.topPercent = topPercent;
	this.maxLength = maxLength;
	this.taxName = getName(id);
	this.filter = f;
	//initlaize with default values for output
	this.coverageLine = taxName+"\tNA\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0";
	ArrayList<Double> pIdents = new ArrayList<Double>();
	ArrayList<Integer> distances = new ArrayList<Integer>();
	ArrayList<String> reads = new ArrayList<String>();
	ArrayList<String> alignments = new ArrayList<String>();
	reads.add("None");
	alignments.add("None");
	pIdents.add(0.0);
	this.readList = reads;
	this.alignmentList = alignments;
	setEditDistanceHistogram(distances);
	setPercentIdentityHistogram(pIdents);
	String s = taxName;
	for(int i = 0;i<=40;i++){
		s+="\t"+0;
	}
	setDamageLine(s);
	
	this.readDistribution = taxName+"\tNA\t0\t0\t0\t0\t0";
	String rldist = taxName;
	for(int i = 25;i<=200;i+=5)
		rldist+="\t0";
	this.readLengthDistribution = rldist;
	readLengthStatistics =taxName+"\t0\t0\t0\t0";
	container.setName(taxName);
	this.additionalEntries =taxName+"\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA";
	this.covPositions = taxName+"\tNA\t0\t0\t0\t0\t0\t0\t0";
	DecimalFormatSymbols otherSymbols = new DecimalFormatSymbols();
	otherSymbols.setDecimalSeparator('.');
	otherSymbols.setGroupingSeparator(','); 
	df =new DecimalFormat("#.###",otherSymbols);
}
//setters
protected void setOriginalNumberOfAlignments(int num){
	this.originalNumberOfAlignments = num;
}
protected void setOriginalNumberOfReads(int num){
	this.originalNumberOfReads = num;
}
protected void setReads(ArrayList<String> list){
	this.readList = list;
}
protected void setAlignments(ArrayList<String> list){
	this.alignmentList = list;
}
protected void setDamageLine(String s){
	this.damageLine = s;
}
// setter for compostion map that retrieves information on read distribution and node composition
protected void processCompositionMap(CompositionMap map){
	if(map.getCompositionMap().keySet().size() >0){//check if some alignments are even left after filtering
		map.calculateStatistics();
		//setReadDist
		String maxReference = getName(map.getMaxID());
		String s = taxName +"\t" + maxReference;;
		for(double d : map.getGenaralStatistics())
			s += "\t" + df.format(d);
		this.readDistribution=s;
		
		//set coverage Line
		HashMap<Integer,Integer> histogram = map.getConverageHistogram();
		String line = taxName + "\t" + maxReference;
		for(int k : histogram.keySet())
			line += "\t" + histogram.get(k);
		this.coverageLine = line;
		this.additionalEntries = taxName+"\t"+map.getTopTenReferences();
		//get coveragePositions
		String covPosLine = taxName+"\t" + maxReference;
		for(String cov : map.getCoveragePositions()){
			covPosLine += "\t"+cov;
		}
		this.covPositions = covPosLine;
		
		map=null; // unassign Map at the end 
	}else{
		this.readDistribution = taxName+"\tNA\t0\t0\t0\t0\t0";
		this.coverageLine = taxName+"\tNA\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0";
		this.covPositions = taxName+"\tNA\t0\t0\t0\t0\t0\t0\t0";
	}
}

protected void setNumberOfReads(int n){
	this.numOfReads=n;
}
// set and process edit distance distribution
protected void setEditDistanceHistogram(ArrayList<Integer> list){
	HashMap<Integer,Integer> histo = new HashMap<Integer,Integer> ();
	for(int i = 0;i < 12;i++){
		histo.put(i, 0);
	}
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
	this.editHistogram = histo;
	list.clear();
	list = null;
}
// process percent identity
protected void setPercentIdentityHistogram(ArrayList<Double> list){
	HashMap<Integer,Integer> histo = new HashMap<Integer,Integer> ();
	for(int i = 0;i < 5; i++){
		histo.put(i, 0);
	}
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
	this.pIdentHistogram =  histo;
	list.clear();
	list = null;
}

protected void setNumMatches(int matches){
	this.numMatches = matches;
}
public void setTurnedOn(boolean b){
	this.turnedOn = b;
}

//getters
public int getTaxID(){
	return taxID;
}
public String getCoveragePositions(){
	return this.covPositions;
}
public boolean wasTurnedOn(){
	return this.turnedOn;
}
public String getFilterLine(){
	String s;
	if(turnedOn)
		s = taxName+"\t"+originalNumberOfReads + "\t" +numOfReads +"\t"+originalNumberOfAlignments+"\t" + numMatches+"\t"+"On";
	else
		s = taxName+"\t"+originalNumberOfReads + "\t" +numOfReads +"\t"+originalNumberOfAlignments+"\t" + numMatches+"\t"+"Off";
	return s;
}
public String getCoverageLine(){
	return this.coverageLine;
}
public String getDamageLine(){
	return this.damageLine;
}
// get Name of class
protected String getName(int taxId){
	String name;
	if(mapReader.getNcbiIdToNameMap().get(taxId) != null)
		name = mapReader.getNcbiIdToNameMap().get(taxId).replace(' ', '_');
	else if(taxId == 0)
		name="NA";
	else
		name = "unassignedName";
	return name;
}
public ArrayList<String> getAlignments(){
	return this.alignmentList;
}
public ArrayList<String> getReads(){
	return this.readList;
}
public String getEditDistanceHistogram(){
	HashMap<Integer,Integer> histo = this.editHistogram;
	return taxName+"\t"+ histo.get(0)+"\t"+histo.get(1)+"\t"+histo.get(2)+"\t"+histo.get(3)+"\t"+histo.get(4)+"\t"+histo.get(5)+"\t"
			+histo.get(6)+"\t"+histo.get(7)+"\t"+histo.get(8)+"\t"+histo.get(9)+"\t"+histo.get(10)+"\t"+histo.get(11);
}
public String getPercentIdentityHistogram(){
	HashMap<Integer,Integer> histo = this.pIdentHistogram;
	return taxName+"\t"+histo.get(0)+"\t"+histo.get(1)+"\t"+histo.get(2)+"\t"+histo.get(3)+"\t"+histo.get(4);
}

protected double getGcContent(String sequence){
	double gcContent = 0;
	char[] chars=sequence.toCharArray();
	for(char c : chars){
		if(c=='g'||c=='G'||c=='c'||c=='C')
			gcContent++;
	}
	if(gcContent !=0){
		gcContent=gcContent/chars.length;
	}
	return gcContent;
}
//calculate read length distribtuion
protected void calculateReadLengthDistribution(){
	DescriptiveStatistics stats = new DescriptiveStatistics();
	HashMap<Integer,Integer> intervals = new HashMap<Integer,Integer>();
	for(int i = 25;i<=200;i+=5)
		intervals.put(i, 0);
	if(lengths.size() == 0 && lengths == null){
		warning.log(Level.WARNING, taxName+" Read Lengths empty");
	}
	else{//check that lengths are set and filled
		boolean warn = false;
		for(int i : lengths){// round to the closest number dividable by 5 to get the intervals
			if(i>=25 && i<=200){//check that lenghts are in interval
				stats.addValue(i);
				int value = intervals.get(round(i, 5));
				value++;
				intervals.replace(round(i, 5), value);
			}else{
				warn = true;
			}
		}
		if(warn) {
			warning.log(Level.WARNING, taxName+"has reads that fall outside the 25-200 target range");
		}
		String rlStat = taxName;
		
		rlStat += "\t"+stats.getMean();
		rlStat += "\t"+stats.getGeometricMean();
		rlStat += "\t"+stats.getPercentile(50);
		rlStat += "\t"+stats.getStandardDeviation();
		
		this.readLengthStatistics = rlStat;
		String rlDist = taxName;
		for(int key:intervals.keySet()){
			rlDist+="\t"+intervals.get(key);
		}
		this.readLengthDistribution = rlDist;
		stats=null;
		lengths.clear();
	}
}
public String getReadLengthDistribution(){
	return this.readLengthDistribution;
}
public String getReadLengthStatistics(){
	return this.readLengthStatistics;
}
public int getNumberOfReads(){
	return this.numOfReads;
}
public HashMap<Integer, HashMap<String, ArrayList<Alignment>>> getTaxonMap(){
	return this.taxonMap;
}
public String getReadDistribution(){
	return this.readDistribution;
}
public String getAdditionalEntries(){
	return this.additionalEntries;
}
private int round(double i, int v){
    return (int) (Math.round(i/v) * v);
}
public void processMatchBlocks(IMatchBlock[] blocks, String name, int length, String sequence){ 

	}// void 
public void process(){ 

}// void 
public void merge(RMA6TaxonProcessor t1){
	lines.addAll(t1.getReads());
	HashMap<Integer, HashMap<String, ArrayList<Alignment>>>merger =  t1.getTaxonMap();
	for(int taxID:merger.keySet()){
		if(!taxonMap.containsKey(taxID)){;
			taxonMap.put(taxID, merger.get(taxID));
		}else{
			HashMap<String,ArrayList<Alignment>> list =  taxonMap.get(taxID);
			HashMap<String,ArrayList<Alignment>> mergedList = merger.get(taxID);
			for(String key: mergedList.keySet()){
				if(!list.containsKey(key)){
					list.put(key, mergedList.get(key));
				}else{
					ArrayList<Alignment> entry = list.get(key);
					entry.addAll(mergedList.get(key));
					list.replace(key,entry);
				}
			}
		}
	}	
}
public void clear(){
	container = null;
	pIdents.clear();
	distances.clear();
}
}// class 
