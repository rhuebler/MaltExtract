package RMA6TaxonProcessor;

import java.io.IOException;
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
 * NonDuplicate Filter for RMA6 Files
 * Due to technical constraints at the moment we only use the highest scoring BLast it when we remove duplicates at the moment
 * Runtime should be longer for duplicate removal due to passing all data one additional time 
 * @author huebler
 *
 */
public class RMA6TaxonNonDuplicateFilter  extends RMA6TaxonProcessor{
	/**
	 * @param int ID, NCBI_MapReader reader, boolean verbose, Logger, log, Logger warning
	 * @return int numMatches, String readDistribution, HashMap EditDistance, HashMap Percent Identity
	 */ 
	public RMA6TaxonNonDuplicateFilter(int id ,double pID, NCBI_MapReader reader, boolean v,Logger log, Logger warning) {
		super(id,pID, reader, v, log, warning);
	}
	private void computeOutput(HashMap<Integer, ArrayList<Alignment>> taxonMap, int taxID){
		HashMap<Integer,Integer> substitutionMap = new HashMap<Integer,Integer>();
		HashMap<Integer,Integer> misMap = new HashMap<Integer,Integer>();
		int numMatches = 0;
		ArrayList<Integer> distances = new ArrayList<Integer>();
		ArrayList<Double> pIdents = new ArrayList<Double>();
		CompositionMap map = new CompositionMap(taxonMap);
		map.process();
		map.markAllDuplicates();
		
		setReadDistribution(map);// first set ReadDistribution on Maximum ID
		taxonMap = map.getCompositionMap();
		int numReads=0;
		for(int key : taxonMap.keySet()){
			for(Alignment entry : taxonMap.get(key)){
				if(!entry.isDuplicate()){
					//get mismatches
					HashMap<Integer, String> alMap = entry.getMismatches();
					if(alMap != null && entry.getMlength() >= 20){//get mismatches per position
						numMatches++;
						for(int l = 0; l< 20; l++){
							if(l < 10){
								if(alMap.containsKey(l))
									if(alMap.get(l).equals("C>T")){
										if(misMap.containsKey(l))
											misMap.replace(l, misMap.get(l)+1);
										else	
											misMap.put(l, 1);	
									}else if(!alMap.get(l).contains("[Nn-]+")){//get all others
										if(substitutionMap.containsKey(l))
											substitutionMap.replace(l, substitutionMap.get(l)+1);
										else
											substitutionMap.put(l, 1);
									}
							}else{
								if(alMap.containsKey(entry.getMlength()+l-20))
									if(alMap.get(entry.getMlength()+l-20).equals("G>A")){
										if(misMap.containsKey(l))
											misMap.replace(l, misMap.get(l)+1);
										else	
											misMap.put(l, 1);
									}else if(!alMap.get(entry.getMlength()+l-20).contains("[Nn-]+")){//get all others
										if(substitutionMap.containsKey(l))
											substitutionMap.replace(l, substitutionMap.get(l)+1);
										else
											substitutionMap.put(l, 1);
									}
								}
							}//for
						}// map != null		
					pIdents.add(entry.getPIdent());
					distances.add(entry.getEditInstance());
					numReads++;
				}
				
			}
			
		}
		setSubstitutionMap(substitutionMap);
		setMisMap(misMap);
		setNumMatches(numMatches);
		setEditDistanceHistogram(distances);
		setPercentIdentityHistogram(pIdents);
		//setSupplementary(supplementary);
		setNumberOfMatches(numReads);
	}	
	@Override
	public void process(String inDir, String fileName, double topPercent, int maxLength){ 
		HashMap<Integer, ArrayList<Alignment>> taxonMap = new HashMap<Integer,ArrayList<Alignment>>();
		// use ReadsIterator to get all Reads assigned to MegantaxID and print top percent to file;
		try(RMA6File rma6File = new RMA6File(inDir+fileName, "r")){
			this.taxName = getName(taxID);
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
					warning.log(Level.WARNING,"TaxName: " +taxName+  " not assigned in File " + fileName+"\n");
				setReadDistribution(new CompositionMap(new HashMap<Integer,ArrayList<Alignment>>()));
				setPercentIdentityHistogram(pIdents);
				setEditDistanceHistogram(distances);
				setNumberOfMatches(0);
				
			}else{
				if(verbose)
					log.log(Level.INFO,"Processing Taxon "+ taxName + " in File " + fileName); 
	
				while(classIt.hasNext()){
				IReadBlock current = classIt.next();
					if(current.getReadLength() <= maxLength || maxLength == 0){
						IMatchBlock block=current.getMatchBlock(0);
						Alignment al = new Alignment();
						al.processText(block.getText().split("\n"));
						al.setReadName(current.getReadName());
						al.setPIdent(block.getPercentIdentity());
						if(minPIdent <= al.getPIdent()){
							if(!taxonMap.containsKey(block.getTaxonId())){
								ArrayList<Alignment> entry =new ArrayList<Alignment>();
								entry.add(al);
								taxonMap.put(block.getTaxonId(), entry);
							}else{
								ArrayList<Alignment> entry = taxonMap.get(block.getTaxonId());
								entry.add(al);
								taxonMap.put(block.getTaxonId(),entry);
							}
						}
				}//if 
			}//while
				classIt.close();
				computeOutput(taxonMap, taxID);
		  }//else		
		}catch(IOException io){
			warning.log(Level.SEVERE,mapReader.getNcbiIdToNameMap().get(taxID), io);
		}
	}	
}
