package RMA6TaxonProcessor;

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

public class SimultaniouslyAncientNon extends RMA6TaxonDamageFilter{

	public SimultaniouslyAncientNon(int id, double pID, NCBI_MapReader reader, boolean v, Logger log, Logger warning) {
		super(id, pID, reader, v, log, warning);
		// TODO Auto-generated constructor stub
	}
	
	public SimultaniouslyAncientNon(int id, double pID, NCBI_MapReader reader, boolean v, Logger log, Logger warning,
			boolean reads) {
		super(id, pID, reader, v, log, warning, reads);
		// TODO Auto-generated constructor stub
	}
	public void process(String inDir, String fileName, double topPercent, int maxLength){ 
		ArrayList<Integer> allDistances = new ArrayList<Integer>();
		ArrayList<Double> allPIdents = new ArrayList<Double>();
		StrainMisMatchContainer allContainer = new StrainMisMatchContainer();
		
		ArrayList<Integer> ancientDistances = new ArrayList<Integer>();
		ArrayList<Double> ancientPIdents = new ArrayList<Double>();
		ArrayList<String> ancientLines = new ArrayList<String>();
		StrainMisMatchContainer ancientContainer = new StrainMisMatchContainer();
		
		this.taxName = getName(taxID);
		int numMatches = 0;
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
			IReadBlockIterator classIt  = new ReadBlockIterator(list, new ReadBlockGetterRMA6(rma6File, true, true, (float) 1.0,(float) 100.00,false,true));
			if(!classIt.hasNext()){ // check if reads are assigned to TaxID if not print to console and skip
				if(verbose)
					warning.log(Level.WARNING,"TaxID: " + taxID +  " not assigned in File " + fileName+"\n");
				setReadDistribution(new CompositionMap(new HashMap<Integer,ArrayList<Alignment>>()));
				setPercentIdentityHistogram(allPIdents);
				setEditDistanceHistogram(allDistances);
				String s = taxName;
				for(int i = 0;i<=40;i++){
					s+="\t"+0;
				}
				setDamageLine(s);
		}else{
			if(verbose)
				log.log(Level.INFO,"Processing Taxon "+taxName+" in File " +fileName); 
			HashMap<Integer, ArrayList<Alignment>> taxonMap = new HashMap<Integer,ArrayList<Alignment>>();
			HashMap<Integer, ArrayList<Alignment>> ancientMap = new HashMap<Integer,ArrayList<Alignment>>();
			int numReads = 0;
			int ancientNumReads = 0;
			while(classIt.hasNext()){
				IReadBlock current = classIt.next();
				boolean higher = false;
				int damage = 0;
				if(current.getReadLength() <= maxLength || maxLength == 0){
					IMatchBlock[] blocks=current.getMatchBlocks();
					int k=0;
					float topScore = current.getMatchBlock(0).getBitScore();
					double pIdent = 0;
					int editDistance=0;
					double ancientPIdent = 0;
					int ancientEditDistance=0;
					for(int i = 0; i< blocks.length;i++){
						if(blocks[i].getBitScore()/topScore < 1-topPercent){
							break;}
						
						Alignment al = new Alignment();
						al.processText(blocks[i].getText().split("\n"));
						al.setPIdent(blocks[i].getPercentIdentity());
						al.setReadName(current.getReadName());
						al.setReadLength(current.getReadLength());
						al.setAcessionNumber(blocks[i].getRefSeqId());	
						if(minPIdent <= al.getPIdent()){ // check for minPercentIdentity
							if(al.getFivePrimeDamage()){
								ancientNumReads++;
								higher = true;
								//get mismatches
								ancientContainer.processAlignment(al);
								if(!ancientMap.containsKey(blocks[i].getTaxonId())){
									ArrayList<Alignment> entry =new ArrayList<Alignment>();
									entry.add(al);
									ancientMap.put(blocks[i].getTaxonId(), entry);
								}else{
									ArrayList<Alignment> entry = ancientMap.get(blocks[i].getTaxonId());
									entry.add(al);
									ancientMap.put(blocks[i].getTaxonId(),entry);
								}
							ancientEditDistance += al.getEditInstance();
							ancientPIdent += al.getPIdent();
							damage++;
							if(wantReads){
								String name = getName(blocks[i].getTaxonId());
								
								ancientLines.add(al.getReadName()+"\t"+"Length:\t"+al.getReadLength()+"\t");
								ancientLines.add(name+"\t"+al.getAccessionNumber()+"\t"+"Start:\t"+al.getStart()+"\t"+"End:\t"+al.getEnd());
								ancientLines.add("Q:\t"+al.getQuery());
								ancientLines.add("A:\t"+al.getAlignment());
								ancientLines.add("R:\t"+al.getReference()+"\n");
							}
							
							}
							//get mismatches
							allContainer.processAlignment(al);
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
						allDistances.add(editDistance/k);
						allPIdents.add(pIdent/k);
					}
					if(higher&&damage !=0){
						ancientNumReads++;
						ancientDistances.add(ancientEditDistance/damage);
						ancientPIdents.add(ancientPIdent/damage);
					}
						
				}// if  
			}// while
				classIt.close();
				//set all non-filter information 
				CompositionMap map = new CompositionMap(taxonMap);
				map.process();
				
				StrainMap strain = new StrainMap(taxName,allContainer,numReads);
				setDamageLine(strain.getLine());
				setNumMatches(numMatches);
				setReadDistribution(map);
				
				setEditDistanceHistogram(allDistances);
				setPercentIdentityHistogram(allPIdents);
				
				//set ancient information
				CompositionMap aMap = new CompositionMap(ancientMap);
				aMap.process();
				
				StrainMap aStrain = new StrainMap(taxName,ancientContainer,ancientNumReads);
				setAncientDamageLine(aStrain.getLine());
				setAncientNumMatches(numMatches);
				setAncientReadDistribution(aMap);
				
				setAncientEditDistanceHistogram(ancientDistances);
				setAncientPercentIdentityHistogram(ancientPIdents);
				setAncientReads(ancientLines);
				
				rma6File.close();
		     }//else
			}catch(Exception e){
				warning.log(Level.SEVERE,mapReader.getNcbiIdToNameMap().get(taxID), e);	
			
			}
		}// void
}
