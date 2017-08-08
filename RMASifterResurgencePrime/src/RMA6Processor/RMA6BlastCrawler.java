package RMA6Processor;

import java.io.IOException;
import java.nio.charset.Charset;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.logging.Level;
import java.util.logging.Logger;

import NCBI_MapReader.NCBI_MapReader;
import NCBI_MapReader.NCBI_TreeReader;
import RMAAlignment.Alignment;
import behaviour.Filter;
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

public class RMA6BlastCrawler {
	/**
	 * Go through whole file and find blast hits that that match the strains of the target species
	 * one line per strain currently only works with Strains that contain species name and won't
	 * work for anything for that is higher than species level
	 */
	//Initilaize attributes
	private String inDir;
	private String fileName;
	private String speciesName;
	private NCBI_MapReader mapReader;
	private String outDir;
	private Logger warning;
	private ArrayList<String> summary= new ArrayList<String>();
	private NCBI_TreeReader treeReader;
	private ArrayList<String> editDistances = new ArrayList<String>();
	private ArrayList<String> percentIdentities = new ArrayList<String>();
	private ArrayList<String> readDistributions = new ArrayList<String>();
	private Filter filter = Filter.CRAWL;
	private ArrayList<String> reads =new ArrayList<String>();
	private ArrayList<String> headers =new ArrayList<String>();
	//set values at construction
	//TODO get Reads
	public RMA6BlastCrawler(String dir, String name, String species, String out, NCBI_MapReader reader ,Logger warning,NCBI_TreeReader treeReader,
			Filter filter){
		this.inDir = dir;
		this.fileName = name;
		this.speciesName = species;
		this.mapReader = reader;
		this.outDir = out;
		this.warning = warning;
		this.treeReader = new NCBI_TreeReader(treeReader);
		this.filter = filter;
	}

	//  Output writers here 
	private void writeReadLengthDistribution(List<String> histo){
		try{
			String header = "Node\tMean\tGeometricMean\tMedian\tStandardDev";
			histo.sort(null);
			histo.add(0,header);
			Path file = Paths.get(outDir+"/crawlResults/readDist/"+fileName+"_"+speciesName.replace(' ', '_')+"_readDist.txt");
			Files.write(file, histo, Charset.forName("UTF-8"));
		}catch(IOException io){
			warning.log(Level.SEVERE,"Cannot write file", io);
		}
	}
	private void writeEditDistance(List<String> histo){
		try{
			String header = "Node\t0\t1\t2\t3\t4\t5\thigher";
			histo.sort(null);
			histo.add(0,header);
			Path file = Paths.get(outDir+"/crawlResults/editDistance/"+fileName+"_"+speciesName.replace(' ', '_')+"_editDistance.txt");
			Files.write(file, histo, Charset.forName("UTF-8"));
		}catch(IOException io){
			warning.log(Level.SEVERE,"Cannot write file", io);
		}
	}
	private void writePercentIdentity(List<String> histo){
		try{
			String header = "Node\t80\t85\t90\t95\t100";
			histo.sort(null);
			histo.add(0,header);
			Path file = Paths.get(outDir+"/crawlResults/percentIdentity/"+fileName+"_"+speciesName.replace(' ', '_')+"_percentIdentity.txt");
			Files.write(file, histo, Charset.forName("UTF-8"));
		}catch(IOException io){
			warning.log(Level.SEVERE,"Cannot write file", io);
		}
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
		Path file = Paths.get(outDir+"/crawlResults/damageMismatch/"+fileName+"_"+speciesName.replace(' ', '_')+"_misMatch.txt");
		Files.write(file, summary, Charset.forName("UTF-8"));
	}catch(IOException io){
		warning.log(Level.SEVERE,"Can't write File" ,io);
	}
 
}
	private void writeReads(List<String>reads,List<String>headers){
		try{
			ArrayList<String> summary = new ArrayList<String>();
			for(int i = 0; i<= reads.size();i++){
				summary.add(headers.get(i));
				summary.add(headers.get(i));
			}
			Path file = Paths.get(outDir+"/crawlResults/reads/"+fileName+"_"+speciesName+"_reads.fa");
			Files.write(file, summary, Charset.forName("UTF-8"));
		}catch(IOException io){
			warning.log(Level.SEVERE,"Cannot write file", io);
		}
	} 	
	// get all keys from a file
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
		    	warning.log(Level.SEVERE,"Can't write File" ,io);io.printStackTrace();
			}
		
		return keys;
	}
	
	
	
	// Process RMA6 file by crawling 
	public void process(){
		// retrieve IDs from file and intilaize StrainMap
		HashMap<Integer, StrainMap> collection = new HashMap<Integer, StrainMap>();
		HashSet<Integer> idsToProcess = new HashSet<Integer>();
		int taxID = mapReader.getNcbiNameToIdMap().get(speciesName);
		idsToProcess.add(taxID);
		for(Integer id : treeReader.getAllStrains(taxID, getAllKeys(inDir+fileName)))
				idsToProcess.add(id);
		idsToProcess.addAll(treeReader.getParents(taxID));
		for(Integer id : idsToProcess){
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
				
				// iterate through all nodes and store information in strain Map to process and retrieve at later use
				while(classIt.hasNext()){
					IReadBlock current = classIt.next();
							IMatchBlock[] blocks = current.getMatchBlocks();
							float topScore = current.getMatchBlock(0).getBitScore();
							for(int i = 0; i< blocks.length;i++){
								if(blocks[i].getBitScore()/topScore < 1-0.01){
									break;}
								if(getName(blocks[i].getTaxonId()).contains(speciesName)){
									if(!reads.contains(current.getReadSequence())&&!headers.contains(current.getReadHeader())){
										reads.add(current.getReadSequence());
										headers.add(current.getReadHeader());
									}
									Alignment al = new Alignment();
									al.setText(blocks[i].getText());
									al.processText();
									al.setPIdent(blocks[i].getPercentIdentity());
									if(collection.containsKey(blocks[i].getTaxonId())){
										StrainMap strain = collection.get(blocks[i].getTaxonId());
										StrainMisMatchContainer container =	strain.getStrainMisMatchContainer();
										container.processAlignment(al);
										strain.setStrainMisMatchContainer(container);
										strain.setNumberOfMatches(strain.getNumberOfMatches()+1);
										collection.replace(blocks[i].getTaxonId(), strain);
							
									}else{
										StrainMisMatchContainer container = new StrainMisMatchContainer(filter);
										container.processAlignment(al);
										StrainMap strain = new StrainMap(getName(blocks[i].getTaxonId()),
										container,1);
										collection.put(blocks[i].getTaxonId(), strain);
											}//else
										}//if
						}// other for 	//matches
				}	
				classIt.close();
				rma6File.close();
			}catch(IOException io){
				warning.log(Level.SEVERE,"Can not locate or read File" ,io);
			}
		}// for all IDs
		for(int key :collection.keySet()){// write output here 
		summary.add(collection.get(key).getLine());
		editDistances.add(collection.get(key).getEditDistanceHistogram());
		percentIdentities.add(collection.get(key).getPercentIdentityHistogram());
		readDistributions.add(collection.get(key).getReadLengthDistribution());
		}
		// write output 
		writeMisMap(summary);
		writeEditDistance(editDistances);
		writePercentIdentity(percentIdentities);
		writeReadLengthDistribution(readDistributions);
		writeReads(reads,headers);
	}
}
