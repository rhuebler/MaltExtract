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
 * Extract all information from one Taxon and save the information in the specified slots to be retrieved
 * in RMA6Processor this taxon processor only processes reads that have at least one match that hints at c-> end substitutions 
 * it will be used to replace the non filtering taxon processor within RMA6Processor when filtering for damaged Reads 
 * @author huebler
 *
 */
public class RMA6TaxonDamageFilter extends RMA6TaxonProcessor{
	protected boolean wantReads = false;
	
	/**
	 * @param int ID, NCBI_MapReader reader, boolean verbose, Logger, log, Logger warning
	 * @return int numMatches, String readDistribution, HashMap EditDistance, HashMap Percent Identity
	 */ 
public RMA6TaxonDamageFilter(int id ,double pID, NCBI_MapReader reader,
		boolean v,Logger log, Logger warning) {
		super(id,pID, reader, v, log, warning);
	}
public RMA6TaxonDamageFilter(int id ,double pID, NCBI_MapReader reader,
		boolean v,Logger log, Logger warning, boolean reads) {
	super(id,pID, reader, v, log, warning);
	this.wantReads =reads;
}
@Override
public void process(String inDir, String fileName, double topPercent, int maxLength){ 
	if(mapReader.getNcbiIdToNameMap().get(taxID) != null)
		this.taxName = getName(taxID);
	ArrayList<Integer> distances = new ArrayList<Integer>();
	ArrayList<Double> pIdents = new ArrayList<Double>();
	ArrayList<String> lines = new ArrayList<String>();
	int numberOfMatches = 0;
	StrainMisMatchContainer container = new StrainMisMatchContainer();
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
		lines.add(taxName);
		IReadBlockIterator classIt  = new ReadBlockIterator(list, new ReadBlockGetterRMA6(rma6File, true, true, (float) 1.0,(float) 100.00,false,true));
		if(!classIt.hasNext()){ // check if reads are assigned to TaxID if not print to console and skip
			if(verbose)
				warning.log(Level.WARNING,"TaxID: " + taxID +  " not assigned in File " + fileName+"\n");
			setReadDistribution(new CompositionMap(new HashMap<Integer,ArrayList<Alignment>>()));
			setPercentIdentityHistogram(pIdents);
			setEditDistanceHistogram(distances);
			setNumberOfReads(0);	
			setReads(lines);
			
			String s = taxName;
			for(int i = 0;i<=40;i++){
				s+="\t"+0;
			}
			setDamageLine(s);
		}else{
			if(verbose)
				log.log(Level.INFO,"Processing Taxon "+ taxName + " in File " + fileName); 
	HashMap<Integer, ArrayList<Alignment>> taxonMap = new HashMap<Integer,ArrayList<Alignment>>();
	int numReads = 0;
		while(classIt.hasNext()){
			IReadBlock current = classIt.next();
			if(current.getReadLength()<= maxLength || maxLength == 0){
				IMatchBlock[] blocks=current.getMatchBlocks();
					int k=0;
				float topScore = blocks[0].getBitScore();
				double pIdent = 0;
				int editDistance = 0;
				int damage=0;
				boolean higher = false;
				for(int i = 0; i< blocks.length;i++){
					if(blocks[i].getBitScore()/topScore <= 1 - topPercent){
						break;}
					Alignment al = new Alignment();
					al.processText(blocks[i].getText().split("\n"));
					al.setPIdent(blocks[i].getPercentIdentity());
					al.setReadName(current.getReadName());
					al.setReadLength(current.getReadLength());
					al.setAcessionNumber(blocks[i].getRefSeqId());
					if(al.getFivePrimeDamage() && minPIdent <= al.getPIdent()){
						numberOfMatches++;
						higher = true;
						pIdent += blocks[i].getPercentIdentity();
						//get mismatches
						container.processAlignment(al);
						if(!taxonMap.containsKey(blocks[i].getTaxonId())){
							ArrayList<Alignment> entry =new ArrayList<Alignment>();
							entry.add(al);
							taxonMap.put(blocks[i].getTaxonId(), entry);
						}else{
							ArrayList<Alignment> entry = taxonMap.get(blocks[i].getTaxonId());
							entry.add(al);
							taxonMap.put(blocks[i].getTaxonId(),entry);
						}
					editDistance += al.getEditInstance();
					pIdent += al.getPIdent();
					damage++;
					k++;
					if(wantReads){
						String name = getName(blocks[i].getTaxonId());
						
						lines.add(al.getReadName()+"\t"+"Length:\t"+al.getReadLength()+"\t");
						lines.add(name+"\t"+al.getAccessionNumber()+"\t"+"Start:\t"+al.getStart()+"\t"+"End:\t"+al.getEnd());
						lines.add("Q:\t"+al.getQuery());
						lines.add("A:\t"+al.getAlignment());
						lines.add("R:\t"+al.getReference()+"\n");
					}
				 }
					
				}
				if(damage !=0 && higher){
					numReads++;
					distances.add(editDistance/k);
					pIdents.add(pIdent/k);
				}	
			}
		}// while
			classIt.close();
			StrainMap strain = new StrainMap(taxName,container,numberOfMatches);
			
			CompositionMap map = new CompositionMap(taxonMap);
			map.process();
			
			setDamageLine(strain.getLine());
			setReadDistribution(map);
			setNumberOfReads(numReads);
			setEditDistanceHistogram(distances);
			setPercentIdentityHistogram(pIdents);
			setReads(lines);
			rma6File.close();
		 }//else
		}catch(IOException io){
			warning.log(Level.SEVERE,mapReader.getNcbiIdToNameMap().get(taxID), io);
		}
	}// void 
}// class 
