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
import java.util.concurrent.ExecutionException;
import java.util.concurrent.Executors;
import java.util.concurrent.ThreadPoolExecutor;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.concurrent.Future;
import NCBI_MapReader.NCBI_MapReader;
import NCBI_MapReader.NCBI_TreeReader;
import RMA6TaxonProcessor.ConcurrentMatchProcessorCrawler;
import RMA6TaxonProcessor.MatchProcessorCrawler;
import RMAAlignment.Alignment;
import behaviour.OutputType;
import jloda.util.ListOfLongs;
import megan.data.IMatchBlock;
import megan.data.IReadBlock;
import megan.data.IReadBlockIterator;
import megan.data.ReadBlockIterator;
import megan.rma6.ClassificationBlockRMA6;
import megan.rma6.RMA6File;
import megan.rma6.ReadBlockGetterRMA6;
import utility.InputParameterProcessor;

public class RMA6BlastCrawler {
	/**
	 * Go through whole file and find blast hits that that match the strains of the target species
	 * one line per strain currently only works with Strains that contain species name and won't
	 * work for anything for that is higher than species level as the name matching will probably fail.
	 * This class is parallelized to some extent to make it a viable option for deeper sequenced data,
	 * however this will still be relatively slow as any read for any alignment on the taxonomic path 
	 * (all strains of that species and all higher nodes up to the root)
	 */
	//Initilaize attributes
	private String inDir;
	private String fileName;
	private String speciesName;
	private NCBI_MapReader mapReader;
	private String outDir;
	private Logger log;
	private Logger warning;
	private ArrayList<String> damageLines= new ArrayList<String>();
	private NCBI_TreeReader treeReader;
	private ArrayList<String> coverageHistograms = new ArrayList<String>();
	private ArrayList<String> coveragePositions = new ArrayList<String>();
	private ArrayList<String> editDistances = new ArrayList<String>();
	private ArrayList<String> percentIdentities = new ArrayList<String>();
	private ArrayList<String> readLengthDistributions = new ArrayList<String>();
	private ArrayList<String> readDistributions = new ArrayList<String>();
	private  HashMap<Integer,ArrayList<Alignment>> concurrentMap  = new HashMap<Integer,ArrayList<Alignment>>();
	private int numThreads;
	private ThreadPoolExecutor executor;
	private boolean singleStranded = false;
	private void destroy(){
		executor.shutdown();
	}
	private  boolean wantReads = false;
	
	private ArrayList<String> reads = new ArrayList<String>();
	//set values at construction
	//TODO get Reads
	public RMA6BlastCrawler(String dir, String name, String species,InputParameterProcessor inProcessor, NCBI_MapReader reader ,Logger log, Logger warning,NCBI_TreeReader treeReader){
		this.inDir = dir;
		this.fileName = name;
		this.speciesName = species;
		this.mapReader = reader;
		this.outDir = inProcessor.getOutDir()+"/crawlResults/";
		this.warning = warning;
		this.treeReader = new NCBI_TreeReader(treeReader);
		this.wantReads = inProcessor.wantReads();
		this.numThreads =inProcessor.getNumThreads();
	}
	private void writeOutput(List<String> histo, String dir,OutputType type,int taxID){
		try{
			String  header ="Taxon";
			String outDir = dir;
			switch (type) {
			case ADDIDTIONALENTRIES: header = "TargetNode\t01\t02\t03\t04\t05\t06\t07\t08\t09\t10";
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
			case EDITDISTANCE: header="Taxon\t0\t1\t2\t3\t4\t5\t6\t7\t8\t9\t10\thigher";
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
		executor=(ThreadPoolExecutor) Executors.newFixedThreadPool(numThreads);
		
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
					String readName = current.getReadName();
					String readSequence = current.getReadSequence();
					int readLength = current.getReadLength();
					float topScore = blocks[0].getBitScore();
					for(int i = 0; i< blocks.length;i++){
						if((blocks[i].getBitScore()/topScore) < 1-0.01){
							break;}	
						String name = getName(blocks[i].getTaxonId());
						
						if(name.contains(speciesName)){
							int cid= blocks[i].getTaxonId();
							//System.out.println("true");
							Alignment al = new Alignment();
							al.setText(blocks[i].getText());
							al.processText();
							al.setTaxID(blocks[i].getTaxonId());
							al.setPIdent(blocks[i].getPercentIdentity());
							al.setReadName(readName);
							al.setReadLength(readLength);
							al.setAcessionNumber(blocks[i].getTextFirstWord());	
							al.setSequence(readSequence);
							al.setTopAlignment(true);
							
							if(concurrentMap.containsKey(cid)) {
								
								ArrayList<Alignment> alist =	concurrentMap.get(cid);
								alist.add(al);
								concurrentMap.replace(cid, alist);
							}else {
								
								ArrayList<Alignment> alist = new ArrayList<Alignment>();
								alist.add(al);
								concurrentMap.put(cid, alist);
							}
							break;
						}
					}
				}	
				classIt.close();
				rma6File.close();
			}catch(IOException io){
				warning.log(Level.SEVERE,"Can not locate or read File" ,io);
			}
		}// for all IDs
		destroy();
		executor=(ThreadPoolExecutor) Executors.newFixedThreadPool(numThreads);
		ArrayList<Future<MatchProcessorCrawler>> futureMPCList = new ArrayList<Future<MatchProcessorCrawler>>(concurrentMap.keySet().size());
		for(int key :concurrentMap.keySet()){
			ConcurrentMatchProcessorCrawler cmpc = new ConcurrentMatchProcessorCrawler(concurrentMap.get(key), key, log, warning, wantReads, mapReader,singleStranded);
			Future<MatchProcessorCrawler> future = executor.submit(cmpc);
			futureMPCList.add(future);
		}
		destroy();
		
		for(Future<MatchProcessorCrawler> future :futureMPCList)
		{// write output here 
			MatchProcessorCrawler matchPC=null;
			try {
				matchPC = future.get();
			} catch (InterruptedException | ExecutionException e) {
				// TODO Auto-generated catch block
				warning.log(Level.SEVERE, "Interuppted Exception", e);
			}
			coverageHistograms.add(matchPC.getCoverageLine());
			coveragePositions.add(matchPC.getCoveragePositions());
			damageLines.add(matchPC.getDamageLine());
			editDistances.add(matchPC.getEditDistanceHistogram());
			percentIdentities.add(matchPC.getPercentIdentityHistogram());
			readLengthDistributions.add(matchPC.getReadLengthDistribution());
			readDistributions.add(matchPC.getReadDistribution());
			reads.addAll(matchPC.getReads());
			matchPC.clear();
		}
		
		writeOutput(coverageHistograms, outDir,OutputType.COVERAGEHISTOGRAM, 0);
		writeOutput(coveragePositions, outDir,OutputType.POS_COVERED, 0);
		writeOutput(editDistances, outDir,OutputType.EDITDISTANCE, 0);
		writeOutput(percentIdentities, outDir,OutputType.PERCENTIDENTITY, 0);
		writeOutput(readDistributions, outDir,OutputType.ALIGNMENTDISTRIBUTION, 0);
		writeOutput(readLengthDistributions, outDir, OutputType.READLENGTH_STATISTICS, 0);
		writeOutput(damageLines, outDir,OutputType.DAMAGE, 0);
		if(wantReads){
			writeOutput(reads,outDir,OutputType.READS,taxID);
		}
	}
}
