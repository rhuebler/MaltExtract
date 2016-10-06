package RMA6TaxonProcessor;

import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.HashMap;

import java.util.logging.Logger;
import NCBI_MapReader.NCBI_MapReader;
import RMAAlignment.CompositionMap;

/**
 * Extract all information from one Taxon and save the information in the specified slots to be retrieved
 * in RMA6Processor is parent class for all filtering Taxonprocessors
 * @author huebler
 *
 */
public class RMA6TaxonProcessor {
	/**
	 * @param int ID, NCBI_MapReader reader, boolean verbose, Logger, log, Logger warning, minPIdent
	 * @return int numMatches, String readDistribution, HashMap EditDistance, HashMap Percent Identity
	 */ 
protected String taxName;	
protected int numOfMatches;
protected String readDistribution;
protected ArrayList<String> supplemantary;
protected NCBI_MapReader mapReader;
protected Integer taxID;
protected double minPIdent;
protected HashMap<Integer,Integer> editHistogram;
protected HashMap<Integer,Integer> pIdentHistogram;
protected boolean verbose;
protected Logger log;
protected Logger warning;
protected ArrayList<String> readList;
protected HashMap<Integer,Integer> misMap;
protected HashMap<Integer,Integer> substitutionMap;
protected int numMatches;
//constructor
public RMA6TaxonProcessor(Integer id, double pID, NCBI_MapReader reader, boolean v, Logger log, Logger warning){
	this.mapReader = reader;
	this.minPIdent = pID;
	this.taxID = id;
	this.verbose = v;
	this.log = log;
	this.warning = warning;
}
//setters
protected void setSubstitutionMap(HashMap<Integer,Integer> map){
	this.substitutionMap = map;
}
protected void setReadDistribution(CompositionMap map){
	DecimalFormat df = new DecimalFormat("#.###");
	if(map != null){
	String maxReference = getName(map.getMaxID());
	String s = taxName +"\t" + maxReference;;
	for(double d : map.getStatistics())
		s += "\t" + df.format(d);
		this.readDistribution=s;
	}else{
		this.readDistribution = mapReader.getNcbiIdToNameMap().get(taxID).replace(' ', '_')+"\tNA\t0\t0\t0\t0\t0\t0\t0\t0";
	}
}

protected void setNumberOfMatches(int n){
	this.numOfMatches=n;
}
protected void setEditDistanceHistogram(ArrayList<Integer> list){
	HashMap<Integer,Integer> histo = new HashMap<Integer,Integer> ();
	for(int i = 0;i < 7;i++){
		histo.put(i, 0);
	}
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
	this.editHistogram = histo;
}
protected void setPercentIdentityHistogram(ArrayList<Double> list){
	HashMap<Integer,Integer> histo = new HashMap<Integer,Integer> ();
	for(int i = 0;i < 5; i++){
		histo.put(i, 0);
	}
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
	this.pIdentHistogram =  histo;
}
protected void setMisMap(HashMap<Integer,Integer>map){
	this.misMap = map;
}
protected void setNumMatches(int matches){
	this.numMatches = matches;
}
//getters
protected String getName(int taxId){
	String name;
	if(mapReader.getNcbiIdToNameMap().get(taxID) != null)
		name = mapReader.getNcbiIdToNameMap().get(taxID).replace(' ', '_');
	else
		name = "unassignedName";
	return name;
}
public HashMap<Integer, Integer>getSubstitutionMap(){
	if( substitutionMap != null){
		int numMatches = this.numMatches;
		
		substitutionMap.put(20, numMatches);
		return this.substitutionMap;
	}else{
		HashMap<Integer,Integer> map = new HashMap<Integer,Integer>();
		map.put(0, 0);
		return map;
	}
	
}
public HashMap<Integer, Integer>getMisMap(){
	if( misMap!=null){
		int numMatches = this.numMatches;
		
		misMap.put(20, numMatches);
		return this.misMap;
	}else{
		HashMap<Integer,Integer> map = new HashMap<Integer,Integer>();
		map.put(0, 0);
		return map;
	}
	
}
public ArrayList<String> getReads(){
	return this.readList;
}
public String getEditDistanceHistogram(){
	HashMap<Integer,Integer> histo = this.editHistogram;
	return taxName+"\t"+ histo.get(0)+"\t"+histo.get(1)+"\t"+histo.get(2)+"\t"+histo.get(3)+"\t"+histo.get(4)+"\t"+histo.get(5)+"\t"+histo.get(6);
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

public int getNumberOfMatches(){
	return this.numOfMatches;
}

public String getReadDistribution(){
	return this.readDistribution;
}

public ArrayList<String> getSupplementary(){
	return this.supplemantary;
}

public void process(String inDir, String fileName, double topPercent, int maxLength){ 

	}// void 
}// class 
