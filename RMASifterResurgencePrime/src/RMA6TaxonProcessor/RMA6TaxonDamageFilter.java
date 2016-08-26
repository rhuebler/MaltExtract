package RMA6TaxonProcessor;

import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.HashMap;

import NCBI_MapReader.NCBI_MapReader;
import RMAAlignment.Alignment;
import RMAAlignment.CompositionMap;
import megan.data.IMatchBlock;
import megan.data.IReadBlock;
import megan.data.IReadBlockIterator;
import megan.rma6.RMA6Connector;
/**
 * Extract all information from one Taxon and save the information in the specified slots to be retrieved
 * in RMA6Processor this taxon processor only processes reads that have at least one match that hints at c-> end substitutions 
 * it will be used to replace the non filtering taxon processor within RMA6Processor when filtering for damaged Reads 
 * @author huebler
 *
 */
public class RMA6TaxonDamageFilter {
private int numOfMatches;
private String readDistribution;
private ArrayList<String> supplemantary;

private void setSupplementary(ArrayList<String> s){
	this.supplemantary=s;
}
private void setReadDistribution(String s){
	this.readDistribution=s;
}
private void setNumberOfMatches(int n){
	this.numOfMatches=n;
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

private double getGcContent(String sequence){
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

public void process(RMA6Connector fileCon, int taxID, String fileName,NCBI_MapReader mapReader, double topPercent, int maxLength){ 
	DecimalFormat df = new DecimalFormat("#.###");
	ArrayList<String> supplemantary = new ArrayList<String>();
	// use ReadsIterator to get all Reads assigned to MegantaxID and print top percent to file;
	try{
	IReadBlockIterator classIt  = fileCon.getReadsIterator("Taxonomy", taxID, (float) 1.0,(float) 100.00,true,true);
	if(!classIt.hasNext()){ // check if reads are assigned to TaxID if not print to console and skip
		System.err.println("TaxID: " + taxID +  " not assigned in File " + fileName+"\n");
		setNumberOfMatches(0);
		}
	// only consider topscoring alignments

	System.out.println("Processing Taxon "+mapReader.getNcbiIdToNameMap().get(taxID)+" in File " + fileName); 
	HashMap<Integer, ArrayList<Alignment>> taxonMap = new HashMap<Integer,ArrayList<Alignment>>();
	int numReads = 0;
		while(classIt.hasNext()){
			IReadBlock current = classIt.next();
			if(current.getReadLength()<= maxLength || maxLength == 0){
				IMatchBlock[] blocks=current.getMatchBlocks();
					int k=0;
				float topScore = blocks[0].getBitScore();
				double pIdent = 0;
				double length = 0; 
				int damage=0;
				for(int i = 0; i< blocks.length;i++){
					if(blocks[i].getBitScore()/topScore <= 1 - topPercent){
						break;}
					
					Alignment al = new Alignment();
					al.processText(blocks[i].getText().split("\n"));
					if(al.getFivePrimeDamage()){
						
						length += al.getMlength();
						pIdent += blocks[i].getPercentIdentity();
						if(!taxonMap.containsKey(blocks[i].getTaxonId())){
							ArrayList<Alignment> entry =new ArrayList<Alignment>();
							entry.add(al);
							taxonMap.put(blocks[i].getTaxonId(), entry);
						}else{
							ArrayList<Alignment> entry = taxonMap.get(blocks[i].getTaxonId());
							entry.add(al);
							taxonMap.put(blocks[i].getTaxonId(),entry);
						}
					
						damage++;
					k++;
					}
				}
				if(damage !=0){
					numReads++;
					supplemantary.add(
						current.getReadName()+"\t"
						+ (int)(length/k)+"\t"
						+ df.format(pIdent/(k)) +"\t"
						+ current.getNumberOfMatches()+"\t"
						+ k +"\t"
						+ damage+'\t'
						+ df.format(getGcContent(current.getReadSequence()))+"\t"
						+ mapReader.getNcbiIdToNameMap().get(taxID).replace(' ', '_'));}
			}// if TODO should I add an else here and what to do 
		}// while
			classIt.close();
			CompositionMap map = new CompositionMap(taxonMap);
			map.process();
			String s = mapReader.getNcbiIdToNameMap().get(taxID).replace(' ', '_') 
						+"\t" + mapReader.getNcbiIdToNameMap().get(map.getMaxID()).replace(' ', '_');
			for(double d : map.getStatistics())
				s += "\t" + df.format(d);
			setReadDistribution(s);
			setNumberOfMatches(numReads);
			setSupplementary(supplemantary);
		}catch(IOException io){
		io.printStackTrace();
		}
	}// void 
}// class 
