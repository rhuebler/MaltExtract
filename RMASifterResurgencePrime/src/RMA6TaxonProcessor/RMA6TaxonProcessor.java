package RMA6TaxonProcessor;

import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.logging.Level;
import java.util.logging.Logger;

import NCBI_MapReader.NCBI_MapReader;
import RMAAlignment.Alignment;
import RMAAlignment.CompositionMap;
import jloda.util.ListOfLongs;
import megan.data.IMatchBlock;
import megan.data.IReadBlock;
import megan.data.IReadBlockIterator;
import megan.data.ReadBlockIterator;
import megan.rma6.ClassificationBlockRMA6;
import megan.rma6.RMA6File;
import megan.rma6.ReadBlockGetterRMA6;
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

protected void setReadDistribution(String s){
	if(s != null){
		this.readDistribution=s;
	}else{
		this.readDistribution = mapReader.getNcbiIdToNameMap().get(taxID).replace(' ', '_')+"\tNA\t0\t0\t0\t0\t0\t0\t0\t0";
	}
}

protected void setNumberOfMatches(int n){
	this.numOfMatches=n;
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
	DecimalFormat df = new DecimalFormat("#.###");
	ArrayList<Integer> distances = new ArrayList<Integer>();
	ArrayList<Double> pIdents = new ArrayList<Double>();
	if(mapReader.getNcbiIdToNameMap().get(taxID) != null)
		this.taxName = mapReader.getNcbiIdToNameMap().get(taxID).replace(' ', '_');
	else
		this.taxName = "unassignedName";
	// use ReadsIterator to get all Reads assigned to MegantaxID and print top percent to file
	try(RMA6File rma6File = new RMA6File(inDir+fileName, "r")){
		ListOfLongs list = new ListOfLongs();
		Long location = rma6File.getFooterSectionRMA6().getStartClassification("Taxonomy");
		if (location != null) {
		   ClassificationBlockRMA6 classificationBlockRMA6 = new ClassificationBlockRMA6("Taxonomy");
		   classificationBlockRMA6.read(location, rma6File.getReader());
		   if (classificationBlockRMA6.getSum(taxID) > 0) {
			   classificationBlockRMA6.readLocations(location, rma6File.getReader(), taxID, list);
		   }
		 }
		HashMap<Integer,Integer> misMap = new HashMap<Integer,Integer>();
		int numMatches = 0;
		IReadBlockIterator classIt  = new ReadBlockIterator(list, new ReadBlockGetterRMA6(rma6File, true, true, (float) 1.0,(float) 100.00,false,true));
		if(!classIt.hasNext()){ // check if reads are assigned to TaxID if not print to console and skip
			if(verbose)
				warning.log(Level.WARNING,"TaxID: " + taxID +  " not assigned in File " + fileName+"\n");
			setReadDistribution(mapReader.getNcbiIdToNameMap().get(taxID).replace(' ', '_')+"\tNA\t0\t0\t0\t0\t0\t0\t0\t0");
			setPercentIdentityHistogram(pIdents);
			setEditDistanceHistogram(distances);
	}else{
		if(verbose)
			log.log(Level.INFO,"Processing Taxon "+mapReader.getNcbiIdToNameMap().get(taxID)+" in File " +fileName); 
		HashMap<Integer, ArrayList<Alignment>> taxonMap = new HashMap<Integer,ArrayList<Alignment>>();
		int numReads = 0;
		while(classIt.hasNext()){
			IReadBlock current = classIt.next();
			boolean higher = false;
			if(current.getReadLength() <= maxLength || maxLength == 0){
				IMatchBlock[] blocks=current.getMatchBlocks();
				int k=0;
				float topScore = current.getMatchBlock(0).getBitScore();
				double pIdent = 0;
				int editDistance=0;
				for(int i = 0; i< blocks.length;i++){
					if(blocks[i].getBitScore()/topScore < 1-topPercent){
						break;}
					
					Alignment al = new Alignment();
					al.processText(blocks[i].getText().split("\n"));
					al.setPIdent(blocks[i].getPercentIdentity());
					//get mismatches
					HashMap<Integer, String> map=al.getMismatches();
					if(map != null && al.getMlength() >= 20){//get mismatches per position
						numMatches++;
						for(int l = 0; l< 20; l++){
							if(l < 10){
								if(map.containsKey(l))
									if(map.get(l).equals("C>T")){
										if(misMap.containsKey(l))
											misMap.replace(l, misMap.get(l)+1);
										else	
											misMap.put(l, 1);	
										}
							}else{
								if(map.containsKey(al.getMlength()+l-20)){
									if(map.get(al.getMlength()+l-20).equals("G>A")){
										if(misMap.containsKey(l))
											misMap.replace(l, misMap.get(l)+1);
										else	
											misMap.put(l, 1);
										}
									}
								}
							}//for
						}// map != null			
					if(minPIdent <= al.getPIdent()){ // check for minPercentIdentity
						higher = true;
						pIdent += al.getPIdent();
						editDistance += al.getEditInstance();
						if(!taxonMap.containsKey(blocks[i].getTaxonId())){
							ArrayList<Alignment> entry =new ArrayList<Alignment>();
							entry.add(al);
							taxonMap.put(blocks[i].getTaxonId(), entry);
						}else{
							ArrayList<Alignment> entry = taxonMap.get(blocks[i].getTaxonId());
							entry.add(al);
							taxonMap.put(blocks[i].getTaxonId(),entry);
						}
						k++;
					}
				}
				if(higher){
					numReads++;
					distances.add(editDistance/k);
					pIdents.add(pIdent/k);
				}
					
			}// if  
		}// while
			classIt.close();
			CompositionMap map = new CompositionMap(taxonMap);
			map.process();
			
			String maxReference;
			if(mapReader.getNcbiIdToNameMap().get(map.getMaxID()) != null)
				maxReference =  mapReader.getNcbiIdToNameMap().get(map.getMaxID()).replace(' ', '_');
			else
				maxReference = "unassinged_reference_name";
			String s = taxName +"\t" + maxReference;;
			for(double d : map.getStatistics())
				s += "\t" + df.format(d);
			this.misMap = misMap;
			this.numMatches = numMatches;
			setReadDistribution(s);
			setNumberOfMatches(numReads);
			setEditDistanceHistogram(distances);
			setPercentIdentityHistogram(pIdents);
			rma6File.close();
	     }//else
		}catch(Exception e){
			warning.log(Level.SEVERE,mapReader.getNcbiIdToNameMap().get(taxID), e);	
		
		}
	}// void 
}// class 
