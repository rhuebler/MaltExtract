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
/**
 * Default filter that just outputs additional proofs without any filtering
 * @author huebler
 *
 */
public class RMA6TaxonNonFilter  extends RMA6TaxonProcessor{
	protected boolean wantReads = false;
	public RMA6TaxonNonFilter(Integer id, double pID, NCBI_MapReader reader, boolean v, Logger log, Logger warning) {
		super(id, pID, reader, v, log, warning);
		// TODO Auto-generated constructor stub
	}
	public RMA6TaxonNonFilter(int id ,double pID, NCBI_MapReader reader,
			boolean v,Logger log, Logger warning, boolean reads) {
		super(id,pID, reader, v, log, warning);
		this.wantReads =reads;
	}
	public void process(String inDir, String fileName, double topPercent, int maxLength){ 
		ArrayList<Integer> distances = new ArrayList<Integer>();
		ArrayList<Double> pIdents = new ArrayList<Double>();
		ArrayList<String> lines = new ArrayList<String>();
		HashMap<Integer,Integer> misMap = new HashMap<Integer,Integer>();
		HashMap<Integer,Integer> substitutionMap = new HashMap<Integer,Integer>();
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
				setPercentIdentityHistogram(pIdents);
				setEditDistanceHistogram(distances);
		}else{
			if(verbose)
				log.log(Level.INFO,"Processing Taxon "+taxName+" in File " +fileName); 
			HashMap<Integer, ArrayList<Alignment>> taxonMap = new HashMap<Integer,ArrayList<Alignment>>();
			int numReads = 0;
			while(classIt.hasNext()){
				IReadBlock current = classIt.next();
				boolean higher = false;
				if(current.getReadLength() <= maxLength || maxLength == 0){
					IMatchBlock[] blocks=current.getMatchBlocks();
					int k=0;
					float topScore = current.getMatchBlock(0).getBitScore();
					double pIdent = 0;
					int editDistance=0;
					for(int i = 0; i< blocks.length;i++){
						if(blocks[i].getBitScore()/topScore < 1-topPercent){
							break;}
						
						Alignment al = new Alignment();
						al.processText(blocks[i].getText().split("\n"));
						al.setPIdent(blocks[i].getPercentIdentity());
							
						if(minPIdent <= al.getPIdent()){ // check for minPercentIdentity
							if(wantReads){
								String name = getName(blocks[i].getTaxonId());
								
								lines.add(al.getReadName()+"\t"+"Length:\t"+al.getReadLength()+"\t");
								lines.add(name+"\t"+al.getAccessionNumber()+"\t"+"Start:\t"+al.getStart()+"\t"+"End:\t"+al.getEnd());
								lines.add("Q:\t"+al.getQuery());
								lines.add("A:\t"+al.getAlignment());
								lines.add("R:\t"+al.getReference()+"\n");
							}
							//get mismatches
							HashMap<Integer, String> map=al.getMismatches();
							if(map != null && al.getMlength() >= 20){//get mismatches per position
								numMatches++;
								for(int l = 0; l< 20; l++){
									if(l < 10){
										if(map.containsKey(l))
											if(map.get(l).equals("C>T")){// get C>T at appropropriate end
												if(misMap.containsKey(l))
													misMap.replace(l, misMap.get(l)+1);
												else	
													misMap.put(l, 1);	
												}else if(!map.get(l).contains("[Nn-]+")){//get all others
													if(substitutionMap.containsKey(l))
														substitutionMap.replace(l, substitutionMap.get(l)+1);
													else
														substitutionMap.put(l, 1);
												}
									}else{
										if(map.containsKey(al.getMlength()+l-20))// get G>A at appropropriate end
											if(map.get(al.getMlength()+l-20).equals("G>A")){
												if(misMap.containsKey(l))
													misMap.replace(l, misMap.get(l)+1);
												else	
													misMap.put(l, 1);
												}else if(!map.get(al.getMlength()+l-20).contains("[Nn-]+")){//get all others
												if(substitutionMap.containsKey(l))
													substitutionMap.replace(l, substitutionMap.get(l)+1);
												else
													substitutionMap.put(l, 1);
											}
										}
									}//for
								}// map != null		
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
						distances.add(editDistance/k);
						pIdents.add(pIdent/k);
					}
						
				}// if  
			}// while
				classIt.close();
				CompositionMap map = new CompositionMap(taxonMap);
				map.process();
				
				setSubstitutionMap(substitutionMap);
				setMisMap(misMap);
				setNumMatches(numMatches);
				setReadDistribution(map);
				setNumberOfMatches(numReads);
				setEditDistanceHistogram(distances);
				setPercentIdentityHistogram(pIdents);
				setReads(lines);
				rma6File.close();
		     }//else
			}catch(Exception e){
				warning.log(Level.SEVERE,mapReader.getNcbiIdToNameMap().get(taxID), e);	
			
			}
		}// void 
}
