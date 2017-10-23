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
import java.util.concurrent.ConcurrentHashMap;
import java.util.logging.Level;
import java.util.logging.Logger;

import NCBI_MapReader.NCBI_MapReader;
import NCBI_MapReader.NCBI_TreeReader;
import RMA6TaxonProcessor.MatchProcessorCrawler;
import behaviour.Filter;
import behaviour.OutputType;
import jloda.util.ListOfLongs;
import megan.data.IMatchBlock;
import megan.data.IReadBlock;
import megan.data.IReadBlockIterator;
import megan.data.ReadBlockIterator;
import megan.rma6.ClassificationBlockRMA6;
import megan.rma6.RMA6File;
import megan.rma6.ReadBlockGetterRMA6;

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
	private ArrayList<String> damageLines= new ArrayList<String>();
	private NCBI_TreeReader treeReader;
	private ArrayList<String> coverageHistograms = new ArrayList<String>();
	private ArrayList<String> coveragePositions = new ArrayList<String>();
	private ArrayList<String> editDistances = new ArrayList<String>();
	private ArrayList<String> percentIdentities = new ArrayList<String>();
	private ArrayList<String> readLengthDistributions = new ArrayList<String>();
	private ArrayList<String> readDistributions = new ArrayList<String>();
	private ConcurrentHashMap<String, MatchProcessorCrawler> concurrentMap = new ConcurrentHashMap<String, MatchProcessorCrawler>();
	private  boolean wantReads = false;
	private Filter filter = Filter.CRAWL;
	private HashMap<String,String> reads = new HashMap<String,String>();
	//set values at construction
	//TODO get Reads
	public RMA6BlastCrawler(String dir, String name, String species, String out, NCBI_MapReader reader ,Logger warning,NCBI_TreeReader treeReader,
			Filter filter, boolean wantReads, int numberOfthreads){
		this.inDir = dir;
		this.fileName = name;
		this.speciesName = species;
		this.mapReader = reader;
		this.outDir = out+"/crawlResults/";
		this.warning = warning;
		this.treeReader = new NCBI_TreeReader(treeReader);
		this.filter = filter;
		this.wantReads = wantReads;
	}
	private void writeOutput(List<String> histo, String dir,OutputType type,int taxID){
		try{
			String  header ="Taxon";
			String outDir = dir;
			switch (type) {
			case ADDIDTIONALENTRIES: header = "TargetNode\t1.0\t0.9\t0.8\t0.7\t0.6\t0.5\t0.4\t0.3\t0.2\t0.1";
									outDir += "readDist/"+fileName+"_additionalNodeEntries"+".txt";
									break;
			case ALIGNMENTS:		outDir += getName(taxID)+".fasta";
									break;
			case ALIGNMENTDISTRIBUTION: header = "Taxon\tReference\tuniquePerReference\tnonStacked\tnonDuplicatesonReference\tTotalAlignmentsOnReference\tReferenceLength";
										outDir+= "readDist/"+fileName+"_alignmentDist"+".txt";
									break;
			case COVERAGEHISTOGRAM: header = "Taxon\tReference\t0\t1\t2\t3\t4\t5\t6\t7\t8\t9\t10\thigher";
									outDir += "coverage/"+fileName+"_coverageHist"+".txt";
									break;
			case DAMAGE: header = "Node";
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
						 	outDir += "damageMismatch/"+fileName+"_damageMismatch"+".txt";
						 			break;
			case EDITDISTANCE: header="Taxon\tReference\t0\t1\t2\t3\t4\t5\t6\t7\t8\t9\t10\thigher";
								outDir += "editDistance/"+fileName+"_editDistance"+".txt";
									break;
			case FILTER: header= "Node\tNumberOfUnfilteredReads\tNumberOfFilteredReads\tNumberOfUnfilteredAlignments\tnumberOfAlignments\tturnedOn?";
						 outDir += "FilterInformation/"+fileName+"_filterTable"+".txt";			
									break;
			case READS:	outDir += "/reads/"+fileName+"_"+getName(taxID)+".fasta";
						break;
			case READLENGTH_DIST: for(int i = 25;i<=200;i+=5){
									header+="\t"+i+"bp";
									}
									outDir+="readDist/"+fileName+"_readLengthDist"+".txt";
									break;
			case READLENGTH_STATISTICS:	header = "Taxon\tMean\tGeometricMean\tMedian\tStandardDev";
									outDir+="readDist/"+fileName+"_readLengthStat"+".txt";
									break;
			
			case PERCENTIDENTITY:	header = "Taxon\t80\t85\t90\t95\t100";
									outDir += "percentIdentity/"+fileName+"_percentIdentity"+".txt";
									break;
			case POS_COVERED:	header = "Taxon\tReference\tAverageCoverage\tCoverge_StandardDeviation\tpercCoveredHigher1\tpercCoveredHigher2\tpercCoveredHigher3\tpercCoveredHigher4\tpercCoveredHigher5";
								outDir += "coverage/"+fileName+"_postionsCovered"+".txt";
									break;
			default:
				warning.log(Level.SEVERE,"No output specified");
				break;
			};
			if(type!=OutputType.READS && type!=OutputType.ALIGNMENTS){
				histo.sort(null);
				histo.add(0,header);
			}
			Path file = Paths.get(outDir);
			Files.write(file, histo, Charset.forName("UTF-8"));
		}catch(IOException io){
			warning.log(Level.SEVERE,"Cannot write file", io);
		}
	}
	//  Output writers here 
	
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
	
	private String getName(int taxId){
		String name;
		if(mapReader.getNcbiIdToNameMap().get(taxId) != null)
			name = mapReader.getNcbiIdToNameMap().get(taxId);
		else
			name = "unassignedName";
		return name;
	}
	
	// Process RMA6 file by crawling 
	public void process(){
		// retrieve IDs from file and intilaize StrainMap
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
//					if(getName(blocks[1].getTaxonId()).contains(speciesName)){
//						if(wantReads){
//							if( !reads.containsKey(current.getReadHeader())){
//								reads.put(current.getReadHeader(),current.getReadSequence());
//							}
//						}	
//					}	
				}	
				classIt.close();
				rma6File.close();
			}catch(IOException io){
				warning.log(Level.SEVERE,"Can not locate or read File" ,io);
			}
		}// for all IDs
//		for(int key :collection.keySet())
		{// write output here 
//			StrainMap map = collection.get(key);
			
//			coverageHistograms.add(map.getCoverageHistogram());
//			coveragePositions.add(map.getCoveragePositions());
//			damageLines.add(collection.get(key).getDamageLine());
//			editDistances.add(collection.get(key).getEditDistanceHistogram());
//			percentIdentities.add(collection.get(key).getPercentIdentityHistogram());
//			readLengthDistributions.add(collection.get(key).getReadLengthDistribution());
//			readDistributions.add(collection.get(key).getReadDistribution());
		}
		ArrayList<String> readsAndHeaders = new ArrayList<String>();
		if(wantReads){
			for(String header: reads.keySet()){
				readsAndHeaders.add(header);
				readsAndHeaders.add(reads.get(header));
				}
		}	
		writeOutput(coverageHistograms, outDir,OutputType.COVERAGEHISTOGRAM, 0);
		writeOutput(coveragePositions, outDir,OutputType.POS_COVERED, 0);
		writeOutput(editDistances, outDir,OutputType.EDITDISTANCE, 0);
		writeOutput(percentIdentities, outDir,OutputType.PERCENTIDENTITY, 0);
		writeOutput(readDistributions, outDir,OutputType.ALIGNMENTDISTRIBUTION, 0);
		writeOutput(readLengthDistributions, outDir, OutputType.READLENGTH_STATISTICS, 0);
		writeOutput(damageLines, outDir,OutputType.DAMAGE, 0);
		if(wantReads){
			writeOutput(readsAndHeaders,outDir,OutputType.READS,taxID);
			}
	}
}
