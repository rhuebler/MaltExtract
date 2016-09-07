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
protected Integer taxID;
protected HashMap<Integer,Integer> editHistogram;
ListOfLongs list;
//constructor
public RMA6TaxonProcessor(Integer id, NCBI_MapReader reader, ListOfLongs l){
	this.mapReader = reader;
	this.taxID = id;
	this.list = l;
}
protected void setEditDistanceHistogram(ArrayList<Integer> list){
	HashMap<Integer,Integer> histo = new HashMap<Integer,Integer> ();
	for(int d : list){
		if(histo.containsKey(d) && d<=5)
		{	int value = histo.get(d);
			histo.put(d, value++);
		}else if(!histo.containsKey(d)&& d<=5){
			histo.put(d, 1);
		}else{
			if(histo.containsKey(6)){
				int value = histo.get(6);
				histo.put(6, value++);
			}else{
				histo.put(6, 1);
			}
		}
		
	}
	this.editHistogram = histo;
}
protected void setPercentIdentityHistogram(ArrayList<Double> list){
	HashMap<Integer,Integer> histo = new HashMap<Integer,Integer> ();
	for(double d: list){
		if(77.5 <= d && d< 82.5){
			if(histo.containsKey(0)){
				int value = histo.get(0);
				histo.put(0, value++);
			}else{
				histo.put(0, 1);
			}
		}else if(82.5 <= d && d< 87.5){
			if(histo.containsKey(1)){
				int value = histo.get(1);
				histo.put(1, value++);
			}else{
				histo.put(1, 1);
			}
		}else if(87.5 <= d && d< 92.5){
			if(histo.containsKey(2)){
				int value = histo.get(2);
				histo.put(2, value++);
			}else{
				histo.put(2, 1);
			}
		}else if(92.5 <= d && d< 97.5){
			if(histo.containsKey(3)){
				int value = histo.get(3);
				histo.put(3, value++);
			}else{
				histo.put(3, 1);
			}
		}else if(97.5 <= d){
			if(histo.containsKey(4)){
				int value = histo.get(4);
				histo.put(4, value++);
			}else{
				histo.put(4, 1);
			}
		}
	}
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

public void process(String inDir, String fileName, double topPercent, int maxLength){ 
	DecimalFormat df = new DecimalFormat("#.###");
	ArrayList<String> supplemantary = new ArrayList<String>();
	ArrayList<Integer> distances = new ArrayList<Integer>();
	ArrayList<Double> pIdents = new ArrayList<Double>();
	// use ReadsIterator to get all Reads assigned to MegantaxID and print top percent to file
	
	try{
	RMA6File rma6File = new RMA6File(inDir+fileName, "r");
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
		System.out.println("Processing Taxon "+mapReader.getNcbiIdToNameMap().get(taxID)+" in File " +fileName); 
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
				int editDistance=0;
				for(int i = 0; i< blocks.length;i++){
					if(blocks[i].getBitScore()/topScore < 1-topPercent){
						break;}
					
					Alignment al = new Alignment();
					al.processText(blocks[i].getText().split("\n"));
					length += al.getMlength();
					al.setPIdent(blocks[i].getPercentIdentity());
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
					if(al.getFivePrimeDamage())// calculate number of damage Reads
						damage++;
					k++;
				}
					numReads++;
					distances.add(editDistance/k);
					pIdents.add(pIdent/k);
					
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
		System.out.println(mapReader.getNcbiIdToNameMap().get(taxID));	
		e.printStackTrace();
		}
	}// void 
}// class 
