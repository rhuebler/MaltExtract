package DatabaseAnalyzer;

import java.io.IOException;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.logging.Level;
import java.util.logging.Logger;

import NCBI_MapReader.NCBI_TreeReader;
import behaviour.Taxas;
import megan.rma6.ClassificationBlockRMA6;
import megan.rma6.RMA6File;
import utility.DataSummaryWriter;

public class ReadDatabaseAnalyzer {
	private String inDir;
	private String fileName;
	private Set<Integer> keySet;
	private Set<Integer> allKeys;
	private Integer readCount;
	private Map<Integer,Integer> assignmentMap;
	private Taxas tax;
	private List<Integer> taxIDs;
	private NCBI_TreeReader  treeReader;
	private String outDir;
	private Logger log;
	private Logger warning;
	private boolean wantMeganSummaries = false;
	// constructor
	public ReadDatabaseAnalyzer(String inDir, String name, Taxas tax, List<Integer> taxIds, NCBI_TreeReader tReader,
			Logger log, Logger warning,String outDir, boolean wantMeganSummaries){
		this.inDir = inDir;
		this.fileName =  name;
		this.outDir = outDir;
		this.taxIDs = taxIds;
		this.tax =  tax;
		this.treeReader = new NCBI_TreeReader(tReader);
		this.log = log;
		this.warning = warning;
		this.wantMeganSummaries = wantMeganSummaries;
		process();
		}
	// getters
	public int getTotalCount(){
		return this.readCount;
	}
	public Set<Integer> getKeySet(){
		return this.keySet;
	}
	
	public Map<Integer,Integer> getAssignmentMap(){
		return this.assignmentMap;
	}
	public String getFileName(){
		return this.fileName;
	}
	// process
	private void process(){
		try{
			log.log(Level.INFO, "Scanning File: "+ fileName);
			if(wantMeganSummaries){
				DataSummaryWriter dsWriter = new DataSummaryWriter(warning);
				dsWriter.writeSummary(inDir, fileName, outDir);
			}
			Map<Integer,Integer> map = new HashMap<Integer,Integer>();
			RMA6File rma6File = new RMA6File(inDir+fileName, "r");
			Long location = rma6File.getFooterSectionRMA6().getStartClassification("Taxonomy");
		    if (location != null) {
		        ClassificationBlockRMA6 cl = new ClassificationBlockRMA6("Taxonomy");
		        cl.read(location, rma6File.getReader());
		        if(tax == Taxas.USER){ //evaluate user taxa list
		        	Set<Integer> idsToProcess = new HashSet<Integer>();
		        	for(Integer taxID : taxIDs){
		        		idsToProcess.add(taxID);
		        		for(Integer id : treeReader.getAllStrains(taxID, cl.getKeySet()))
		        			if(!taxIDs.contains(id))
		        				idsToProcess.add(id);
		        	}
		        	this.keySet =  idsToProcess;
		        	this.allKeys = cl.getKeySet();
		        }else{//if no taxa list provided use all taxa 
		        	this.keySet = cl.getKeySet();
		        	this.allKeys = cl.getKeySet();
				}
		        this.readCount = (int) rma6File.getFooterSectionRMA6().getNumberOfReads();// read in all counts
		        for(int key : keySet){
		        		if(allKeys.contains(key))
		        			map.put(key, cl.getSum(key));
		        		else
		        			map.put(key, 0);
		        }
			this.assignmentMap = map;
		    }else{
		    	warning.log(Level.SEVERE,fileName+" has no taxonomy block");
		    	map.put(0, 0);
		    	this.assignmentMap = map;
		    }
			rma6File.close();
		}catch(IOException io){
			warning.log(Level.SEVERE,"Cannot open File",io);
		}
	}
}
