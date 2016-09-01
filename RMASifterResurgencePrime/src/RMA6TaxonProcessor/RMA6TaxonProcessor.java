package RMA6TaxonProcessor;


import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;

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
protected int numOfMatches;
protected String readDistribution;
protected ArrayList<String> supplemantary;
protected NCBI_MapReader mapReader;
protected int taxID;

public RMA6TaxonProcessor(int id, NCBI_MapReader reader){
	this.mapReader = reader;
	this.taxID = id;
}

protected void setSupplementary(ArrayList<String> s){
	if(s != null){
		this.supplemantary = s;
	}else{
		ArrayList<String> list = new ArrayList<String>();
		list.add("0\t0\t0\t0\t0\t0\t"+mapReader.getNcbiIdToNameMap().get(taxID).replace(' ', '_'));
		this.supplemantary =list;
	}
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

public void process(RMA6File rma6File, String fileName, double topPercent, int maxLength){ 
	DecimalFormat df = new DecimalFormat("#.###");
	ArrayList<String> supplemantary = new ArrayList<String>();
	// use ReadsIterator to get all Reads assigned to MegantaxID and print top percent to file
	
	try{
	ClassificationBlockRMA6 block = new ClassificationBlockRMA6("Taxonomy");
    long start = rma6File.getFooterSectionRMA6().getStartClassification("Taxonomy");
    block.read(start, rma6File.getReader());
    ListOfLongs list = new ListOfLongs();
        if (block.getSum(taxID) > 0) {
            block.readLocations(start, rma6File.getReader(), taxID, list);
        }
	IReadBlockIterator classIt  = new ReadBlockIterator(list, new ReadBlockGetterRMA6(rma6File, true, true, (float) 1.0,(float) 100.00,false,true));
	if(!classIt.hasNext()){ // check if reads are assigned to TaxID if not print to console and skip
		//System.err.println("TaxID: " + taxID +  " not assigned in File " + fileName+"\n");
		setReadDistribution(mapReader.getNcbiIdToNameMap().get(taxID).replace(' ', '_')+"\tNA\t0\t0\t0\t0\t0\t0\t0\t0");
		setSupplementary(new ArrayList<String>(Arrays.asList("0\t0\t0\t0\t0\t0\t"+mapReader.getNcbiIdToNameMap().get(taxID).replace(' ', '_'))));// in case of taxID not being supported add empty Line
		setSupplementary(supplemantary);
	}else{
		String taxName;
		if(mapReader.getNcbiIdToNameMap().get(taxID) != null)
			taxName = mapReader.getNcbiIdToNameMap().get(taxID).replace(' ', '_');
		else
			taxName = "unassignedName";
		//System.out.println("Processing Taxon "+mapReader.getNcbiIdToNameMap().get(taxID)+" in File " +fileName); 
		HashMap<Integer, ArrayList<Alignment>> taxonMap = new HashMap<Integer,ArrayList<Alignment>>();
		int numReads = 0;
		while(classIt.hasNext()){
			IReadBlock current = classIt.next();
			if(current.getReadLength() <= maxLength || maxLength == 0){
				IMatchBlock[] blocks=current.getMatchBlocks();
				int k=0;
				float topScore = current.getMatchBlock(0).getBitScore();
				double pIdent = 0;
				double length = 0; 
				int damage=0;
				for(int i = 0; i< blocks.length;i++){
					if(blocks[i].getBitScore()/topScore < 1-topPercent){
						break;}
					
					pIdent += blocks[i].getPercentIdentity();
					Alignment al = new Alignment();
					al.processText(blocks[i].getText().split("\n"));
					length += al.getMlength();
					if(!taxonMap.containsKey(blocks[i].getTaxonId())){
						ArrayList<Alignment> entry =new ArrayList<Alignment>();
						entry.add(al);
						taxonMap.put(blocks[i].getTaxonId(), entry);
					}else{
						ArrayList<Alignment> entry = taxonMap.get(blocks[i].getTaxonId());
						entry.add(al);
						taxonMap.put(blocks[i].getTaxonId(),entry);
					}
					if(al.getFivePrimeDamage())// calculate number of damage Reads
						damage++;
					k++;
				}
					numReads++;
					supplemantary.add(
						current.getReadName()+"\t"
						+ (int)(length/k)+"\t"
						+ df.format(pIdent/(k)) +"\t"
						+ current.getNumberOfMatches()+"\t"
						+ k +"\t"
						+ damage+'\t'
						+ df.format(getGcContent(current.getReadSequence()))+"\t"
						+ taxName);
			}// if TODO should I add an else here and what to do with it 
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
				s+="\t" + df.format(d);
			setReadDistribution(s);
			setNumberOfMatches(numReads);
			setSupplementary(supplemantary);
			rma6File.close();
	     }//else
		}catch(Exception e){
		e.printStackTrace();
		}
	}// void 
}// class 
