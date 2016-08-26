package RMA6TaxonProcessor;

import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;

import NCBI_MapReader.NCBI_MapReader;
import RMAAlignment.Alignment;
import RMAAlignment.CompositionMap;
import megan.data.IMatchBlock;
import megan.data.IReadBlock;
import megan.data.IReadBlockIterator;
import megan.rma6.RMA6Connector;
/**
 * NonDuplicate Filter for RMA6 Files
 * Due to technical constraints at the moment we only use the highest scoring BLast it when we remove duplicates at the moment
 * Runtime should be longer for duplicate removal due to passing all data one additional time 
 * @author huebler
 *
 */
public class RMA6TaxonNonDuplicateFilter {
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
		ArrayList<String> supplementary = new ArrayList<String>();
		HashMap<Integer, ArrayList<Alignment>> taxonMap = new HashMap<Integer,ArrayList<Alignment>>();
		// use ReadsIterator to get all Reads assigned to MegantaxID and print top percent to file;
		try{
			IReadBlockIterator classIt  = fileCon.getReadsIterator("Taxonomy", taxID, (float) 1.0,(float) 100.00,true,true);
			if(!classIt.hasNext()){ // check if reads are assigned to TaxID if not print to console and skip could potentially only happen if some genus is unavailable 
				System.err.println("TaxID: " + taxID +  " not assigned in File " + fileName+"\n");
				setReadDistribution(mapReader.getNcbiIdToNameMap().get(taxID).replace(' ', '_')+"\tNA\t0\t0\t0\t0\t0\t0\t0\t0");
				setSupplementary(new ArrayList<String>(Arrays.asList("0\t0\t0\t0\t0\t0\t"+mapReader.getNcbiIdToNameMap().get(taxID).replace(' ', '_'))));// in case of taxID not being supported add empty Line
				setNumberOfMatches(0);
				}
			System.out.println("Processing Taxon "+mapReader.getNcbiIdToNameMap().get(taxID)+" in File " + fileName); 
	
			while(classIt.hasNext()){
				IReadBlock current = classIt.next();
				if(current.getReadLength() <= maxLength || maxLength == 0){
					IMatchBlock block=current.getMatchBlock(0);
						Alignment al = new Alignment();
						al.processText(block.getText().split("\n"));
						al.setReadName(current.getReadName());
						al.setPIdent(block.getPercentIdentity());
						if(!taxonMap.containsKey(block.getTaxonId())){
							ArrayList<Alignment> entry =new ArrayList<Alignment>();
							entry.add(al);
							taxonMap.put(block.getTaxonId(), entry);
						}else{
							ArrayList<Alignment> entry = taxonMap.get(block.getTaxonId());
							entry.add(al);
							taxonMap.put(block.getTaxonId(),entry);
						}
				}//if 
			}//while
		}catch(IOException io){
			io.printStackTrace();
		}
		CompositionMap map = new CompositionMap(taxonMap);
		map.process();
		map.markAllDuplicates();
		// first set ReadDistribution on Maximum ID
		String s = mapReader.getNcbiIdToNameMap().get(taxID).replace(' ', '_') + "\t" 
					+ mapReader.getNcbiIdToNameMap().get(map.getMaxID()).replace(' ', '_');
		for(double d : map.getStatistics())
			s+="\t" + df.format(d);
		setReadDistribution(s);
		taxonMap = map.getCompositionMap();
		int numReads=0;
		for(int key : taxonMap.keySet()){
			for(Alignment entry : taxonMap.get(key)){
				if(!entry.isDuplicate()){
					int damage=0;
					if(entry.getFivePrimeDamage())
						damage=1;
					supplementary.add(entry.getReadName()+"\t"
								+ entry.getPIdent()+"\t"
								+ df.format(entry.getNumGaps()) +"\t"
								+ 1 + "\t"
								+ 1 + "\t"
								+ damage + "\t"
								+ df.format(getGcContent(entry.getQuery()))+"\t"
								+ mapReader.getNcbiIdToNameMap().get(taxID).replace(' ', '_'));
					numReads++;
				}
				
			}
			
		}
		setSupplementary(supplementary);
		setNumberOfMatches(numReads);
	}	
}
