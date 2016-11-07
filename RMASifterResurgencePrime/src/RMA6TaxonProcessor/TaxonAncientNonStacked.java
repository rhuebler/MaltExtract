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
import strainMap.StrainMap;
import strainMap.StrainMisMatchContainer;
/**
 * NonDuplicate Filter for RMA6 Files
 * Due to technical constraints at the moment we only use the highest scoring BLast it when we remove duplicates at the moment
 * Runtime should be longer for duplicate removal due to passing all data one additional time 
 * @author huebler
 *
 */
public class TaxonAncientNonStacked  extends RMA6TaxonDamageFilter{
	/**
	 * @param int ID, NCBI_MapReader reader, boolean verbose, Logger, log, Logger warning
	 * @return int numMatches, String readDistribution, HashMap EditDistance, HashMap Percent Identity
	 */ 
	public TaxonAncientNonStacked(int id ,double pID, NCBI_MapReader reader,
			boolean v,Logger log, Logger warning) {
		super(id,pID, reader, v, log, warning);
	}
	public TaxonAncientNonStacked(int id ,double pID, NCBI_MapReader reader,
			boolean v,Logger log, Logger warning, boolean reads) {
		super(id,pID, reader, v, log, warning,reads);
	}
	
	private void computeOutput(HashMap<Integer, ArrayList<Alignment>> taxonMap, int taxID){
		ArrayList<Integer> distances = new ArrayList<Integer>();
		ArrayList<Double> pIdents = new ArrayList<Double>();
		ArrayList<String> lines = new ArrayList<String>();
		StrainMisMatchContainer container = new StrainMisMatchContainer();
		int numMatches = 0;
		lines.add(taxName);
		CompositionMap map = new CompositionMap(taxonMap);
		map.process();
		map.markAllDuplicates();
		setReadDistribution(map);
		taxonMap = map.getCompositionMap();
		int numReads=0;
		for(int key : taxonMap.keySet()){
			for(Alignment entry : taxonMap.get(key)){
				if(!entry.isDuplicate()){
					String name;
					if(mapReader.getNcbiIdToNameMap().get(key) != null)
						name =  mapReader.getNcbiIdToNameMap().get(key).replace(' ', '_');
					else
						name = "unassinged_reference_name";
					lines.add(entry.getReadName()+"\t"+"Length:\t"+entry.getReadLength()+"\t");
					lines.add(name+"\t"+entry.getAccessionNumber()+"\t"+"Start:\t"+entry.getStart()+"\t"+"End:\t"+entry.getEnd());
					lines.add("Q:\t"+entry.getQuery());
					lines.add("A:\t"+entry.getAlignment());
					lines.add("R:\t"+entry.getReference()+"\n");
					//get mismatches
					numMatches++;
					container.processAlignment(entry);
					pIdents.add(entry.getPIdent());
					distances.add(entry.getEditInstance());
					numReads++;
				}
				
			}
			
		}
		StrainMap strain = new StrainMap(taxName,container,numMatches);
		setDamageLine(strain.getLine());
		setNumberOfReads(numReads);
		setEditDistanceHistogram(distances);
		setPercentIdentityHistogram(pIdents);
		
	}	
	@Override
	public void process(String inDir, String fileName, double topPercent, int maxLength){ 
		
			this.taxName = getName(taxID);
		
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
				ArrayList<String> lines = new ArrayList<String>();
				lines.add(taxName);
				if(verbose)
					warning.log(Level.WARNING,"TaxID: " + taxID +  " not assigned in File " + fileName+"\n");
				setReadDistribution(new CompositionMap(new HashMap<Integer,ArrayList<Alignment>>()));
				setPercentIdentityHistogram(pIdents);
				setEditDistanceHistogram(distances);
				String s = taxName;
				for(int i = 0;i<=40;i++){
					s+="\t"+0;
				}
				setDamageLine(s);
				setReads(lines);
			}else{
				if(verbose)
					log.log(Level.INFO,"Processing Taxon "+ taxName + " in File " + fileName); 
			
				while(classIt.hasNext()){
					IReadBlock current = classIt.next();
					if(current.getReadLength() <= maxLength || maxLength == 0){
						IMatchBlock block = current.getMatchBlock(0);
							Alignment al = new Alignment();
							al.processText(block.getText().split("\n"));
							al.setReadName(current.getReadName());
							al.setPIdent(block.getPercentIdentity());
							if(al.getFivePrimeDamage() && minPIdent <= al.getPIdent()){
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
			warning.log(Level.SEVERE,mapReader.getNcbiIdToNameMap().get(taxID), io);
		}
	}	
}
