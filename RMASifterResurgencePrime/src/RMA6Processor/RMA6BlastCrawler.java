package RMA6Processor;

import java.io.IOException;
import java.nio.charset.Charset;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Set;

import NCBI_MapReader.NCBI_MapReader;
import RMAAlignment.Alignment;
import jloda.util.ListOfLongs;
import megan.data.IMatchBlock;
import megan.data.IReadBlock;
import megan.data.IReadBlockIterator;
import megan.data.ReadBlockIterator;
import megan.rma6.ClassificationBlockRMA6;
import megan.rma6.RMA6File;
import megan.rma6.ReadBlockGetterRMA6;
import strainMap.StrainMap;

public class RMA6BlastCrawler {
	/**
	 * Go through whole file and find blast hits that that match the strains of the target species
	 */
	
	//TODO what am I missing here 
	private String inDir;
	private String fileName;
	private String speciesName;
	private NCBI_MapReader mapReader;
	private String outDir;
	private static ArrayList<String> summary= new ArrayList<String>();
	public RMA6BlastCrawler(String dir, String name, String species, String out, NCBI_MapReader reader ){
		this.inDir = dir;
		this.fileName = name;
		this.speciesName = species;
		this.mapReader = reader;
		this.outDir = out;
	}
	private void prepareMisMap(HashMap<Integer,Integer> misMap,HashMap<Integer,Integer> substitutionMap, String name, int numberOfMatches){
		String part1 = name;
		String part2 = "";
		for(int i = 0;i < 20; i++){
			if(misMap.containsKey(i)){
				part1 += "\t" + misMap.get(i);
			}else{	
				part1 += "\t" + 0;	
			}
			if(substitutionMap.containsKey(i)){
				part2 += "\t" + substitutionMap.get(i);
			}else{
				part2 += "\t" + 0;
			}	
		}
		part1 += part2 + "\t" + numberOfMatches;
		summary.add(part1);
	}
	private String getName(int taxId){
		String name;
		if(mapReader.getNcbiIdToNameMap().get(taxId) != null)
			name = mapReader.getNcbiIdToNameMap().get(taxId);
		else
			name = "unassignedName";
		return name;
	}
	private void writeMisMap(ArrayList<String> summary){
		String header = "Node";
		String header_part2 ="";
		for(int i = 0; i < 20; i++){
			if(i<10){
				header+="\t"+"C>T_"+(i+1);
				header_part2+="\t"+"D>V(11Substitution)_"+(i+1);
			}else{
				header+="\t"+"G>A_"+(i+1);
				header_part2+="\t"+"H>B(11Substitution)_"+(i+1);
			}
		}
		header+=header_part2+"\tconsidered_Matches";
		summary.sort(null);
		summary.add(0,header);
	try{
		Path file = Paths.get(outDir+fileName+speciesName.replace(' ', '_')+"_misMatch.txt");
		Files.write(file, summary, Charset.forName("UTF-8"));
	}catch(IOException io){
		io.printStackTrace();
	}
 
}
	private Set<Integer> getAllKeys(String fileName){
		Set<Integer> keys = null;
		try(RMA6File rma6File = new RMA6File(fileName, "r")){
			Long location = rma6File.getFooterSectionRMA6().getStartClassification("Taxonomy");
		    if (location != null) {
		        ClassificationBlockRMA6 classificationBlockRMA6 = new ClassificationBlockRMA6("Taxonomy");
		        classificationBlockRMA6.read(location, rma6File.getReader());
		        keys = classificationBlockRMA6.getKeySet();// get all assigned IDs in a file 
		    }
		    rma6File.close();
		    }catch(IOException io){
				io.printStackTrace();
			}
		
		return keys;
	}
	public void process(){
		HashMap<Integer, StrainMap> collection = new HashMap<Integer, StrainMap>();
		mapReader = new NCBI_MapReader("/Users/huebler/git/RMA6Sifter/RMASifterResurgencePrime/resources/");
		for(Integer id : getAllKeys(inDir+fileName)){
			try(RMA6File rma6File = new RMA6File(inDir+fileName, "r")){
				ListOfLongs list = new ListOfLongs();
				Long location = rma6File.getFooterSectionRMA6().getStartClassification("Taxonomy");
				if (location != null) {
				   ClassificationBlockRMA6 classificationBlockRMA6 = new ClassificationBlockRMA6("Taxonomy");
				   classificationBlockRMA6.read(location, rma6File.getReader());
				   if (classificationBlockRMA6.getSum(id) > 0) {
					   classificationBlockRMA6.readLocations(location, rma6File.getReader(), id, list);
				   }
				 }
				IReadBlockIterator classIt  = new ReadBlockIterator(list, new ReadBlockGetterRMA6(rma6File, true, true, (float) 1.0,(float) 100.00,false,true));
				// get all stuff
				while(classIt.hasNext()){
					IReadBlock current = classIt.next();
							IMatchBlock[] blocks = current.getMatchBlocks();
							float topScore = current.getMatchBlock(0).getBitScore();
							for(int i = 0; i< blocks.length;i++){
								if(blocks[i].getBitScore()/topScore < 1-0.01){
									break;}
								if(getName(blocks[i].getTaxonId()).contains(speciesName)){
				
									Alignment al = new Alignment();
									al.processText(blocks[i].getText().split("\n"));
									HashMap<Integer, String> map=al.getMismatches();
									if(collection.containsKey(blocks[i].getTaxonId())){
										StrainMap strain = collection.get(blocks[i].getTaxonId());
										HashMap<Integer,Integer> misMap = strain.getMisMap();
										HashMap<Integer,Integer> substitutionMap = strain.getSubstitutionMap();
										if(map != null && al.getMlength() >= 20){//get mismatches per position
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
														}//else
													}//for
												}//if Map != Null
										strain.setMisMap(misMap);
										strain.setSubstitutionMap(substitutionMap);
										strain.setNumberOfMatches(strain.getNumberOfMatches()+1);
										collection.replace(blocks[i].getTaxonId(), strain);
									}else{
										HashMap<Integer,Integer> misMap = new HashMap<Integer,Integer>();
										HashMap<Integer,Integer> substitutionMap = new HashMap<Integer,Integer>();
										if(map != null && al.getMlength() >= 20){//get mismatches per position
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
												}//if Map != Null
											StrainMap strain = new StrainMap(mapReader.getNcbiIdToNameMap().get(blocks[i].getTaxonId()),
													misMap, substitutionMap,1);
											collection.put(blocks[i].getTaxonId(), strain);
											}//else
										}//if
							}// other for 	//matches
				}	
				classIt.close();
				rma6File.close();
			}catch(IOException io){
				io.printStackTrace();
			}
		}// for all IDs
		for(int key :collection.keySet())
			prepareMisMap(collection.get(key).getMisMap(),collection.get(key).getSubstitutionMap(),
					collection.get(key).getName(),collection.get(key).getNumberOfMatches());
		writeMisMap(summary);
	}
}
