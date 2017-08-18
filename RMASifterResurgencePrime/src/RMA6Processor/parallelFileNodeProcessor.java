package RMA6Processor;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.ThreadPoolExecutor;
import java.util.logging.Level;
import java.util.logging.Logger;

import NCBI_MapReader.NCBI_MapReader;
import NCBI_MapReader.NCBI_TreeReader;
import RMA6OutputProcessor.RMA6OutputProcessor;
import RMA6TaxonProcessor.ConcurrentNodeProcessor;
import RMA6TaxonProcessor.NodeProcessor;
import behaviour.Filter;
import behaviour.Taxas;
import megan.rma6.ClassificationBlockRMA6;
import megan.rma6.RMA6File;
import utility.DataSummaryWriter;
/**
 * take care of extracting all the information for one RMA6 file and a List of taxons
 * contains functions to write its own output for supplementary information
 * uses several taxons processors to process RMAFiles in different ways to reflect 
 * the user filter choice and make adding new filters more accessible 
 * @author huebler
 *
 */
public class parallelFileNodeProcessor {
	/**
	 * @param String indDir, String fileName, String outDir, NCBI_MapReader, NCBI_TreeReader, Set<Integer>, Filter behave,
	 * Taxas taxas, boolean verbose, int maxLength, Logger log, Logger warning
	 * @return HashMap<Integer,Integer> HashMap<taxID,numOfMatches>
	 * @throws none thrown all caught
	 */
	//set attributes
	private HashMap<String, HashMap<Integer,Integer>> overallSum = new HashMap<String, HashMap<Integer,Integer>>();
	private HashMap<String, HashMap<Integer,Integer>> ancientSum = new HashMap<String, HashMap<Integer,Integer>>();
	
	private String outDir;
	private String[] inFiles;
	private NCBI_MapReader mapReader;
	private NCBI_TreeReader treeReader;
	private Set<Integer> containedIDs;
	private Filter behave;
	private int maxLength;
	private double minPIdent;
	private Taxas taxas;
	private boolean verbose;
	private HashMap<String,Integer> totalCounts = new HashMap<String,Integer>();
	private Logger log;
	private Logger warning;
	private boolean reads;
	private boolean alignments;
	private ThreadPoolExecutor executor;
	private int threads = 1;
	private double minComplexity;
	private boolean wantMeganSummaries;
	private boolean turnOffDestacking;
	private boolean turnOffDeDuping;
	
	// constructor and intilaize attributes
	public parallelFileNodeProcessor(String[] inFiles, String outDir, NCBI_MapReader mapReader,
			NCBI_TreeReader treeReader, int maxLength, double pIdent, Filter b, Taxas t, boolean verbose,
			Logger log, Logger warning, boolean readInf, double minCompl, boolean alignment, boolean wantMeganSummaries, boolean turnOffDestacking, boolean dedupOff) {
		this.inFiles = inFiles;
		this.outDir = outDir;
		this.mapReader = mapReader;
		this.treeReader = new NCBI_TreeReader(treeReader);
		this.behave = b;
		this.maxLength = maxLength;
		this.minPIdent = pIdent;
		this.taxas = t;
		this.verbose = verbose;
		this.log = log;
		this.warning = warning;
		this.alignments = alignment;
		this.reads = readInf;
		this.minComplexity = minCompl;
		this.wantMeganSummaries = wantMeganSummaries;
		this.turnOffDestacking = turnOffDestacking;
		this.turnOffDeDuping = dedupOff;
	}
	

	// private utility functions
	private void destroy(){
		executor.shutdown();
	}
	
	//getters
	private Set<Integer> getAllKeys(String inDir, String fileName){
		Set<Integer> keys = null;
		try(RMA6File rma6File = new RMA6File(inDir+fileName, "r")){
			Long location = rma6File.getFooterSectionRMA6().getStartClassification("Taxonomy");
		    if (location != null) {
		        ClassificationBlockRMA6 classificationBlockRMA6 = new ClassificationBlockRMA6("Taxonomy");
		        classificationBlockRMA6.read(location, rma6File.getReader());
		        keys = classificationBlockRMA6.getKeySet();// get all assigned IDs in a file 
		    }
		   totalCounts.put(fileName,(int) rma6File.getFooterSectionRMA6().getNumberOfReads());
		    rma6File.close();
		    }catch(IOException io){
				io.printStackTrace();
			}
		
		return keys;
	}
public void process(List<Integer>taxIDs, double topPercent) {// processing through file 
	executor=(ThreadPoolExecutor) Executors.newFixedThreadPool(threads);
	ArrayList<Future<NodeProcessor>> results =  new ArrayList<Future<NodeProcessor>>();
	for(String file: inFiles){
		 File f = new File(file);
		 String inDir = f.getParent() +"/";
		 String fileName = f.getName();
		log.log(Level.INFO,"Reading File: " +inDir+ fileName );
		if(wantMeganSummaries){
			DataSummaryWriter dsWriter = new DataSummaryWriter(warning);
			dsWriter.writeSummary(inDir, fileName, outDir);
		}
		Set<Integer> keys = getAllKeys(f.getParent() +"/",f.getName());
		Set<Integer> idsToProcess = new HashSet<Integer>();
		if(taxas == Taxas.USER){// use user specified taxas 
			for(Integer taxID : taxIDs){
				idsToProcess.add(taxID);
				idsToProcess.addAll(treeReader.getAllStrains(taxID, keys));
			}
		}else if(taxas == Taxas.ALL){// or use all taxas here has to be triggered at program initialization
				idsToProcess.addAll(keys);
		}
		containedIDs.addAll(idsToProcess);// write down contained IDs
		for(Integer id : idsToProcess){
			NodeProcessor nodeProcessor = new NodeProcessor(id, minPIdent, mapReader, verbose,log, warning,reads,behave, minComplexity,alignments,turnOffDestacking,turnOffDeDuping);
			ConcurrentNodeProcessor task = new ConcurrentNodeProcessor(nodeProcessor,inDir, fileName, topPercent, maxLength);
			Future<NodeProcessor> future = executor.submit(task);
			results.add(future);
		}//TaxIDs	
	  }
	destroy();
	HashMap<String,HashMap<Integer,NodeProcessor>> sortedNodes = new HashMap<String,HashMap<Integer,NodeProcessor>>();
	for(Future<NodeProcessor> futureNode:results){
		try {
			NodeProcessor node = futureNode.get();
			String fileName = node.getFileName();
			if(sortedNodes.containsKey(fileName)){
				HashMap<Integer,NodeProcessor> map = sortedNodes.get(fileName);
				map.put(node.getTaxId(),node);
				sortedNodes.replace(fileName, map);
			}else{
				HashMap<Integer,NodeProcessor> map = new HashMap<Integer,NodeProcessor>();
				map.put(node.getTaxId(),node);
				sortedNodes.replace(fileName, map);
			}
		} catch (InterruptedException | ExecutionException e) {
			
			warning.log(Level.SEVERE,"Severe Error",e);
			
		}
	  }
	}
	public void calculateRunOutpu(HashMap<String,HashMap<Integer,NodeProcessor>> sortedNodes){
		for(String fileName:sortedNodes.keySet()){
			RMA6OutputProcessor outProcessor = new RMA6OutputProcessor(fileName, outDir,mapReader,warning, behave,alignments, reads);// initilaize output processor here 
			outProcessor.process(sortedNodes.get(fileName));
				if(behave == Filter.NON_ANCIENT){
					overallSum.put(fileName, outProcessor.getSumLine());
					ancientSum.put(fileName, outProcessor.getAncientLine());
				}else if(behave == Filter.ANCIENT){
					ancientSum.put(fileName, outProcessor.getAncientLine());
				}else if(behave == Filter.NON){
					overallSum.put(fileName, outProcessor.getSumLine());
				}
			} 
	}
 }
