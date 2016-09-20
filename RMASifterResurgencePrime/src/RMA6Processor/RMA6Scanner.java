package RMA6Processor;

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

/**
 * Class returns the keyset of a RMA6File
 * which means it returns the names of all Nodes that have an assigned Read
 * it returns the everything as a hash map with the taxId serving as a key 
 * and the number of assigned Reads as result
 * @author huebler
 *
 */
public class RMA6Scanner {
	/**
	 * @param String inDir, String FileName, String outDir, List<Integer> taxIDs
	 * Taxas tax,TreeReader reader, Logger log, Logger warning
	 * @return ap<Integer,Integer> assignmentMap
	 * @throws none thrown all caught
	 */
	// attributes
	private String inDir;
	private String fileName;
	private Set<Integer> keySet;
	private Set<Integer> allKeys;
	private Integer readCount;
	private Map<Integer,Integer> assignmentMap;
	private Taxas tax;
	private List<Integer> taxIDs;
	private NCBI_TreeReader  treeReader;
	private Logger log;
	private Logger warning;
	// constructor
	public RMA6Scanner(String inDir, String name, Taxas tax, List<Integer> taxIds, NCBI_TreeReader tReader,
			Logger log, Logger warning){
		this.inDir = inDir;
		this.fileName =  name;
		this.taxIDs = taxIds;
		this.tax =  tax;
		this.treeReader = new NCBI_TreeReader(tReader);
		this.log = log;
		this.warning = warning;
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
			Map<Integer,Integer> map = new HashMap<Integer,Integer>();
			RMA6File rma6File = new RMA6File(inDir+fileName, "r");
			Long location = rma6File.getFooterSectionRMA6().getStartClassification("Taxonomy");
		    if (location != null) {
		        ClassificationBlockRMA6 cl = new ClassificationBlockRMA6("Taxonomy");
		        cl.read(location, rma6File.getReader());
		        if(tax == Taxas.USER){
		        	Set<Integer> idsToProcess = new HashSet<Integer>();
		        	for(Integer taxID : taxIDs){
		        		idsToProcess.add(taxID);
		        		for(Integer id : treeReader.getStrains(taxID, cl.getKeySet()))
		        			if(!taxIDs.contains(id))
		        				idsToProcess.add(id);
		        	}
		        	this.keySet =  idsToProcess;
		        	this.allKeys = cl.getKeySet();
		        }else{
		        	this.keySet = cl.getKeySet();
		        	this.allKeys = cl.getKeySet();
				}
		        int sum = 0;
		        for(int key : allKeys){
		        	if(keySet.contains(key))
		        		map.put(key, cl.getSum(key));
		        	sum += cl.getSum(key);
		        }
			this.assignmentMap = map;
			this.readCount =  sum;
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
