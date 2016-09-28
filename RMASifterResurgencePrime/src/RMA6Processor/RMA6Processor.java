package RMA6Processor;

import java.io.File;
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
import RMA6TaxonProcessor.RMA6TaxonDamageFilter;
import RMA6TaxonProcessor.RMA6TaxonNonDuplicateFilter;
import RMA6TaxonProcessor.RMA6TaxonProcessor;
import RMA6TaxonProcessor.TaxonAncientNonStacked;
import behaviour.Filter;
import behaviour.Taxas;
import megan.rma6.ClassificationBlockRMA6;
import megan.rma6.RMA6File;
/**
 * take care of extracting all the information for one RMA6 file and a List of taxons
 * contains functions to write its own output for supplementary information
 * uses several taxons processors to process RMAFiles in different ways to reflect 
 * the user filter choice and make adding new filters more accessible 
 * @author huebler
 *
 */
public class RMA6Processor {
	/**
	 * @param String indDir, String fileName, String outDir, NCBI_MapReader, NCBI_TreeReader, Set<Integer>, Filter behave,
	 * Taxas taxas, boolean verbose, int maxLength, Logger log, Logger warning
	 * @return HashMap<Integer,Integer> HashMap<taxID,numOfMatches>
	 * @throws none thrown all caught
	 */
	private HashMap<Integer,Integer> overallSum;
	private ArrayList<String> readDist;
	private ArrayList<String> supplement;
	private String outDir;
	private String fileName;
	private String inDir;
	private NCBI_MapReader mapReader;
	private NCBI_TreeReader treeReader;
	private Set<Integer> containedIDs;
	private Filter behave;
	private int maxLength;
	private double minPIdent;
	private Taxas taxas;
	private boolean verbose;
	private int totalCount;
	private Logger log;
	private Logger warning;
	private boolean reads;
	// constructor
	public RMA6Processor(String inDir, String fileName, String outDir, NCBI_MapReader mapReader,
			NCBI_TreeReader treeReader, int maxLength, double pIdent, Filter b, Taxas t, boolean verbose,
			Logger log, Logger warning) {
		this.inDir = inDir;
		this.outDir = outDir;
		this.fileName = fileName;
		this.mapReader = mapReader;
		this.treeReader = new NCBI_TreeReader(treeReader);
		this.behave = b;
		this.maxLength = maxLength;
		this.minPIdent = pIdent;
		this.taxas = t;
		this.verbose = verbose;
		this.log = log;
		this.warning = warning;
	}
	public RMA6Processor(String inDir, String fileName, String outDir, NCBI_MapReader mapReader,
			NCBI_TreeReader treeReader, int maxLength, double pIdent, Filter b, Taxas t, boolean verbose,
			Logger log, Logger warning, boolean readInf) {
		this.inDir = inDir;
		this.outDir = outDir;
		this.fileName = fileName;
		this.mapReader = mapReader;
		this.treeReader = new NCBI_TreeReader(treeReader);
		this.behave = b;
		this.maxLength = maxLength;
		this.minPIdent = pIdent;
		this.taxas = t;
		this.verbose = verbose;
		this.log = log;
		this.warning = warning;
		this.reads = readInf;
	}
	
	//setters
	private Set<Integer> getAllKeys(){
		Set<Integer> keys = null;
		try(RMA6File rma6File = new RMA6File(inDir+fileName, "r")){
			Long location = rma6File.getFooterSectionRMA6().getStartClassification("Taxonomy");
		    if (location != null) {
		        ClassificationBlockRMA6 classificationBlockRMA6 = new ClassificationBlockRMA6("Taxonomy");
		        classificationBlockRMA6.read(location, rma6File.getReader());
		        keys = classificationBlockRMA6.getKeySet();// get all assigned IDs in a file 
		    }
		    this.totalCount = (int) rma6File.getFooterSectionRMA6().getNumberOfReads();
		    rma6File.close();
		    }catch(IOException io){
				io.printStackTrace();
			}
		
		return keys;
	}
	
	private void setContainedIDs(Set<Integer> set){
		this.containedIDs = set;
	}
	private void setSumLine(HashMap<Integer,Integer> list)
	{	this.overallSum = list;
	}
	// private utility functions
	private void writeBlastHits(ArrayList<String> summary, int taxID){
		try{
			if(summary.size()>1){
			String name;
			if(mapReader.getNcbiIdToNameMap().get(taxID) != null)
				name = mapReader.getNcbiIdToNameMap().get(taxID).replace(' ', '_');
			else
				name = "unassingned name";
			Path file = Paths.get(outDir+"reads/"+fileName.substring(0,fileName.length()-4)+"/"+name+".txt");
			Files.write(file, summary, Charset.forName("UTF-8"));
			}
		}catch(IOException io){
			warning.log(Level.SEVERE,"Cannot write file", io);
		}
	}
	private void writeReadDist(List<String> summary, String fileName){
		try{
			summary.sort(null);
			String header = "Taxon\tReference\tMeanReadDistance\tMedianReadDistance\tVarianceReadDistance\tStandardDeviationReadDistance\tuniquePerReference\tnonDuplicatesonReference\tTotalReadsOnReference\tReferenceLength";
			summary.add(0, header);
			Path file = Paths.get(outDir+"/readDist/"+fileName+"_readDist"+".txt");
			Files.write(file, summary, Charset.forName("UTF-8"));
		}catch(IOException io){
			warning.log(Level.SEVERE,"Cannot write file", io);
		}
		this.readDist = null; //delete data to save space after were done potentially the gc should take of it
	}
	private void writeEditDistance(List<String> histo){
		try{
			String header = "Node\t0\t1\t2\t3\t4\t5\thigher";
			histo.sort(null);
			histo.add(0,header);
			Path file = Paths.get(outDir+"/editDistance/"+fileName+"_editDistance"+".txt");
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
			Path file = Paths.get(outDir+"/percentIdentity/"+fileName+"_percentIdentity"+".txt");
			Files.write(file, histo, Charset.forName("UTF-8"));
		}catch(IOException io){
			warning.log(Level.SEVERE,"Cannot write file", io);
		}
	} 
	
	//getter
	public int getTotalCount(){
		return this.totalCount;
		
	}
	public ArrayList<String> getReadDistribution(){
		return this.readDist;
	}

	public ArrayList<String> getSupplementary(){
		return this.supplement;
	}
	public HashMap<Integer,Integer> getSumLine(){
		return this.overallSum;
	}
	public Set<Integer> getContainedIDs(){
		return this.containedIDs;
	}
	public String getfileName(){
		return this.fileName;
	}
// processing 
public void process(List<Integer>taxIDs, double topPercent) {
	log.log(Level.INFO,"Reading File: " +inDir+fileName);
	HashMap<Integer,Integer> overallSum = new HashMap<Integer,Integer>();
	ArrayList<String> editDistance = new ArrayList<String>();
	ArrayList<String> percentIdentity = new ArrayList<String>();
	ArrayList<String> readDistribution = new ArrayList<String>();
	if(reads)
		new File(outDir+"reads/"+fileName.substring(0, fileName.length()-4)+"/").mkdirs();
	Set<Integer> keys = getAllKeys();
	Set<Integer> idsToProcess = new HashSet<Integer>();
   // treeReader here to avoid synchronization issues 
	if(taxas == Taxas.USER){
		for(Integer taxID : taxIDs){
			idsToProcess.add(taxID);
			for(Integer id : treeReader.getStrains(taxID, keys))
				if(!taxIDs.contains(id))
					idsToProcess.add(id);
		}
	}
	else if(taxas == Taxas.ALL){
		idsToProcess.addAll(keys);
	}
	setContainedIDs(idsToProcess);
		for(Integer id : idsToProcess){
			RMA6TaxonProcessor taxProcessor = null;
			if(behave == Filter.NON){// change to all
				taxProcessor = new RMA6TaxonProcessor(id, minPIdent, mapReader, verbose,log, warning);
			}else if(behave == Filter.ANCIENT){
				 if(reads)
					 taxProcessor = new RMA6TaxonDamageFilter(id, minPIdent, mapReader, verbose,log, warning, reads);
				 else	 
					taxProcessor = new RMA6TaxonDamageFilter(id, minPIdent, mapReader, verbose,log, warning);
			}else if(behave == Filter.NONDUPLICATES){
				 taxProcessor = new RMA6TaxonNonDuplicateFilter(id, minPIdent, mapReader, verbose, log, warning);
			}else if(behave == Filter.ALL){
				if(reads)
					taxProcessor = new TaxonAncientNonStacked(id, minPIdent, mapReader, verbose, log, warning, reads);
				else
					taxProcessor = new TaxonAncientNonStacked(id, minPIdent, mapReader, verbose, log, warning);
			}
			ConcurrentRMA6TaxonProcessor task = new ConcurrentRMA6TaxonProcessor(taxProcessor,inDir, fileName, topPercent, maxLength);
			Future<RMA6TaxonProcessor> future = executor.submit(task);
			results.put(id, future);
			taxProcessor.process(inDir, fileName, topPercent, maxLength);
			overallSum.put(id,taxProcessor.getNumberOfMatches());
			readDistribution.add(taxProcessor.getReadDistribution());
			editDistance.add(taxProcessor.getEditDistanceHistogram());
			percentIdentity.add(taxProcessor.getPercentIdentityHistogram());
			if((behave == Filter.ALL && reads )|| (behave == Filter.ANCIENT && reads)){
				writeBlastHits(taxProcessor.getReads(),id);
			}
	  }//TaxIDs	
	setSumLine(overallSum); // set number of assigned Reads to overall file summary
	writeReadDist(readDistribution,fileName); // RMA6Processor now saves its own output 
	writeEditDistance(editDistance);
	writePercentIdentity(percentIdentity);
    }
 }
