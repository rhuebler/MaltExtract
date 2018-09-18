package RMA6Processor;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.logging.Level;
import java.util.logging.Logger;
import NCBI_MapReader.NCBI_MapReader;
import NCBI_MapReader.NCBI_TreeReader;
import RMA6OutputProcessor.RMA6OutputProcessor;
import RMA6TaxonProcessor.NodeProcessor;
import behaviour.Filter;
import behaviour.Taxas;
import megan.rma6.ClassificationBlockRMA6;
import megan.rma6.RMA6File;
import utility.DataSummaryWriter;
import utility.InputParameterProcessor;
/**
 * Take care of extracting all the information for one RMA6 file and a List of Taxon IDS
 * uses RMA6OutputProcessor to write output and use NodeProcessor to process each node. 
 * Also stores some information for Summary Writer. One RMA6Processor can run independently per thread.
 * For that NCBI Map Reader and NCBI Tree Reader are copied at initialization
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
	private HashMap<Integer,Integer> ancientSum;
	private String outDir;
	private String fileName;
	private String inDir;
	private NCBI_MapReader mapReader;
	private NCBI_TreeReader treeReader;
	private Set<Integer> containedIDs = new HashSet<Integer>();
	private Filter behave;
	private int maxLength;
	private double minPIdent;
	private Taxas taxas;
	private boolean verbose;
	private int totalCount;
	private Logger log;
	private Logger warning;
	private boolean reads;
	private double minComplexity;
	private boolean alignments;
	private boolean wantMeganSummaries;
	private boolean turnOffDestacking;
	private boolean turnOffDeDuping;
	private List<Integer>taxIDs; 
	private double topPercent;
	private boolean downsample;
	private boolean useAllAlignments;
	private boolean singleStranded;
	// constructor
	public RMA6Processor(InputParameterProcessor inputParameterProcessor, String inDir, String fileName, ArrayList<Integer> taxIDs, NCBI_MapReader mapReader,NCBI_TreeReader treeReader,
			Logger log, Logger warning) {
		this.inDir = inDir;
		this.fileName = fileName;
		this.taxIDs = taxIDs;
		this.log = log;
		this.warning = warning;
		this.mapReader = mapReader;
		this.outDir = inputParameterProcessor.getOutDir();
		this.treeReader = new NCBI_TreeReader(treeReader);
		this.behave = inputParameterProcessor.getFilter();
		this.maxLength = inputParameterProcessor.getMaxLength();
		this.minPIdent = inputParameterProcessor.getMinPIdent();
		this.taxas = inputParameterProcessor.getTaxas();
		this.verbose = inputParameterProcessor.isVerbose();
		this.minComplexity = inputParameterProcessor.getMinComplexity();
		this.reads = inputParameterProcessor.wantReads();
		this.wantMeganSummaries = inputParameterProcessor.wantMeganSummaries();
		this.topPercent = inputParameterProcessor.getTopPercent();
		this.alignments = inputParameterProcessor.getBlastHits();
		this.turnOffDestacking =  inputParameterProcessor.turnDestackingOff();
		this.turnOffDeDuping = inputParameterProcessor.getDeDupOff();
		this.downsample = inputParameterProcessor.downsampling();
		this.useAllAlignments = inputParameterProcessor.useAllAlignments();
		this.singleStranded = inputParameterProcessor.isSingleStranded();
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
		    	warning.log(Level.SEVERE, "Could not retrieve all keys for file: "+fileName, io);
			}
		
		return keys;
	}
	
	private void setContainedIDs(Set<Integer> set){
		if(set!=null && set.size() != 0) {
			containedIDs.addAll(set);
		}else {
			containedIDs.add(-1);
			totalCount = 0;
		}
	}
	private void setSumLine(HashMap<Integer,Integer> list){	
		this.overallSum = list;
	}
	// private utility functions
	
	//getter
	public HashMap<Integer,Integer> getAncientLine(){
		return this.ancientSum;
	}
	
	public int getTotalCount(){
		return this.totalCount;
		
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

public void process() {// processing 
		log.log(Level.INFO,"Reading File: " +inDir+fileName);
		if(wantMeganSummaries){
			DataSummaryWriter dsWriter = new DataSummaryWriter(warning);
			dsWriter.writeSummary(inDir, fileName, outDir);
		}
		Set<Integer> keys = getAllKeys();
		Set<Integer> idsToProcess = new HashSet<Integer>();
	   // treeReader here to avoid synchronization issues 
		if(taxas == Taxas.USER && behave!= Filter.SRNA){
			for(Integer taxID : taxIDs){
				idsToProcess.add(taxID);
				idsToProcess.addAll(treeReader.getAllStrains(taxID, keys));
			}
		}else if(taxas == Taxas.USER && behave== Filter.SRNA) {
			idsToProcess.addAll(treeReader.getLowestContainedIds(keys));
		}else if(taxas == Taxas.ALL){
			idsToProcess.addAll(keys);
		}
		setContainedIDs(idsToProcess);
		HashMap<Integer,NodeProcessor> results =  new HashMap<Integer,NodeProcessor>();
			for(Integer id : idsToProcess){
				NodeProcessor nodeProcessor = new NodeProcessor(id, minPIdent, mapReader, verbose,log, warning,reads,behave, minComplexity,alignments,turnOffDestacking,turnOffDeDuping,downsample, useAllAlignments, singleStranded);
				nodeProcessor.process(inDir, fileName, topPercent, maxLength);
				results.put(id, nodeProcessor);
		  }//TaxIDs	
	
		RMA6OutputProcessor outProcessor = new RMA6OutputProcessor(fileName, outDir,mapReader,warning, behave, alignments, reads);
		outProcessor.process(results);
		switch(behave) {
			default:
				setSumLine(outProcessor.getSumLine());
				ancientSum = outProcessor.getAncientLine();
				break;
			case NON:
				setSumLine(outProcessor.getSumLine());
				break;
			case ANCIENT:
				ancientSum = outProcessor.getAncientLine();
				break;
		}
		log.log(Level.INFO, "Finished File: "+fileName);
    }
 }
