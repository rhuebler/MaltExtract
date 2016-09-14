package RMA6TaxonProcessor;

import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
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
public class TaxonAncientNonStacked  extends RMA6TaxonProcessor{
	
	public TaxonAncientNonStacked(int id, NCBI_MapReader reader, boolean v) {
		super(id, reader, v);
	}
	
	private void computeOutput(HashMap<Integer, ArrayList<Alignment>> taxonMap, int taxID){
		ArrayList<Integer> distances = new ArrayList<Integer>();
		ArrayList<Double> pIdents = new ArrayList<Double>();
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
					pIdents.add(entry.getPIdent());
					distances.add(entry.getEditInstance());
					numReads++;
				}
				
			}
			
		}
		setEditDistanceHistogram(distances);
		setPercentIdentityHistogram(pIdents);
		//setSupplementary(supplementary);
		setNumberOfMatches(numReads);
	}	
	@Override
	public void process(String inDir, String fileName, double topPercent, int maxLength){ 
		if(mapReader.getNcbiIdToNameMap().get(taxID) != null)
			taxName = mapReader.getNcbiIdToNameMap().get(taxID).replace(' ', '_');
		else
			taxName = "unassingned name";
		HashMap<Integer, ArrayList<Alignment>> taxonMap = new HashMap<Integer,ArrayList<Alignment>>();
		// use ReadsIterator to get all Reads assigned to MegantaxID and print top percent to file;
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
			IReadBlockIterator classIt  = new ReadBlockIterator(list, new ReadBlockGetterRMA6(rma6File, true, true, (float) 1.0,(float) 100.00,false,true));
			if(!classIt.hasNext()){ // check if reads are assigned to TaxID if not print to console and skip could potentially only happen if some genus is unavailable 
				ArrayList<Integer> distances = new ArrayList<Integer>();
				ArrayList<Double> pIdents = new ArrayList<Double>();
				if(verbose)
					System.err.println("TaxID: " + taxID +  " not assigned in File " + fileName+"\n");
				setReadDistribution(taxName + "\tNA\t0\t0\t0\t0\t0\t0\t0\t0");
				setPercentIdentityHistogram(pIdents);
				setEditDistanceHistogram(distances);
				setNumberOfMatches(0);
				
			}else{
				if(verbose){
					System.out.println("Processing Taxon "+ taxName + " in File " + fileName); 
				}
				while(classIt.hasNext()){
					IReadBlock current = classIt.next();
					if(current.getReadLength() <= maxLength || maxLength == 0){
						IMatchBlock block = current.getMatchBlock(0);
							Alignment al = new Alignment();
							al.processText(block.getText().split("\n"));
							al.setReadName(current.getReadName());
							al.setPIdent(block.getPercentIdentity());
							if(al.getFivePrimeDamage()){
								if(!taxonMap.containsKey(block.getTaxonId())){
									ArrayList<Alignment> entry =new ArrayList<Alignment>();
									entry.add(al);
									taxonMap.put(block.getTaxonId(), entry);
								}else{
									ArrayList<Alignment> entry = taxonMap.get(block.getTaxonId());
									entry.add(al);
									taxonMap.put(block.getTaxonId(),entry);
									}
							}//5' damage	
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