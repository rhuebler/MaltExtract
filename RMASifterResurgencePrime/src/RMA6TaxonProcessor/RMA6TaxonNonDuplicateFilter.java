package RMA6TaxonProcessor;

import java.io.IOException;
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
 * NonDuplicate Filter for RMA6 Files
 * Due to technical constraints at the moment we only use the highest scoring BLast it when we remove duplicates at the moment
 * Runtime should be longer for duplicate removal due to passing all data one additional time 
 * @author huebler
 *
 */
public class RMA6TaxonNonDuplicateFilter  extends RMA6TaxonProcessor{
	
	public RMA6TaxonNonDuplicateFilter(int id, NCBI_MapReader reader) {
		super(id, reader);
		// TODO Auto-generated constructor stub
	}
	private NCBI_MapReader mapReader;
	
	private void computeOutput(HashMap<Integer, ArrayList<Alignment>> taxonMap, int taxID){
		ArrayList<String> supplementary = new ArrayList<String>();
		DecimalFormat df = new DecimalFormat("#.###");
		String taxName;
		if(mapReader.getNcbiIdToNameMap().get(taxID) != null)
			taxName = mapReader.getNcbiIdToNameMap().get(taxID).replace(' ', '_');
		else
			taxName = "unassingned name";
		
		CompositionMap map = new CompositionMap(taxonMap);
		map.process();
		map.markAllDuplicates();
		// first set ReadDistribution on Maximum ID
		String maxReference;
		if(mapReader.getNcbiIdToNameMap().get(map.getMaxID()) != null)
			maxReference =  mapReader.getNcbiIdToNameMap().get(map.getMaxID()).replace(' ', '_');
		else
			maxReference = "unassinged_reference_name";
		String s = taxName + "\t" 
					+ maxReference;
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
								+ taxName);
					numReads++;
				}
				
			}
			
		}
		setSupplementary(supplementary);
		setNumberOfMatches(numReads);
	}	
	@Override
	public void process(RMA6File rma6File, String fileName, double topPercent, int maxLength){ 
		HashMap<Integer, ArrayList<Alignment>> taxonMap = new HashMap<Integer,ArrayList<Alignment>>();
		// use ReadsIterator to get all Reads assigned to MegantaxID and print top percent to file;
		try{
		final ClassificationBlockRMA6 cBlock = new ClassificationBlockRMA6("Taxonomy");
	    final long start = rma6File.getFooterSectionRMA6().getStartClassification("Taxonomy");
	    cBlock.read(start, rma6File.getReader());
	    final ListOfLongs list = new ListOfLongs();
	        if (cBlock.getSum(taxID) > 0) {
	        	cBlock.readLocations(start, rma6File.getReader(), taxID, list);
	        }
			IReadBlockIterator classIt  = new ReadBlockIterator(list, new ReadBlockGetterRMA6(rma6File, true, true, (float) 1.0,(float) 100.00,false,true));
			if(!classIt.hasNext()){ // check if reads are assigned to TaxID if not print to console and skip could potentially only happen if some genus is unavailable 
				System.err.println("TaxID: " + taxID +  " not assigned in File " + fileName+"\n");
				setReadDistribution(mapReader.getNcbiIdToNameMap().get(taxID).replace(' ', '_')+"\tNA\t0\t0\t0\t0\t0\t0\t0\t0");
				setSupplementary(new ArrayList<String>(Arrays.asList("0\t0\t0\t0\t0\t0\t"+mapReader.getNcbiIdToNameMap().get(taxID).replace(' ', '_'))));// in case of taxID not being supported add empty Line
				setNumberOfMatches(0);
				
			}else{
				//System.out.println("Processing Taxon "+mapReader.getNcbiIdToNameMap().get(taxID)+" in File " + fileName); 
	
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
				classIt.close();
				computeOutput(taxonMap, taxID);
		  }//else		
		}catch(IOException io){
			io.printStackTrace();
		}
	}	
}
