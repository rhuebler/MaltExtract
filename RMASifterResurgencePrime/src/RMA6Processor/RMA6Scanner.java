package RMA6Processor;

import java.io.IOException;
import java.util.HashMap;
import java.util.Map;
import java.util.Set;

import megan.rma6.RMA6Connector;

/**
 * Class returns the keyset of a RMA6File
 * which means it returns the names of all Nodes that have an assigned Read
 * it returns the everything as a hash map with the taxId serving as a key 
 * and the number of assigned Reads as result
 * @author huebler
 *
 */
public class RMA6Scanner {
	// attributes
	private String inDir;
	private String fileName;
	private Set<Integer> keySet;
	private Map<Integer,Integer> assignmentMap;
	// constructor
	public RMA6Scanner(String inDir, String name){
		this.inDir = inDir;
		this.fileName =  name;
		process();
		}
	// getters
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
		try{System.out.println("Scanning File: "+ fileName);
			Map<Integer,Integer> map = new HashMap<Integer,Integer>();
			RMA6Connector fileCon = new RMA6Connector(inDir+fileName);
			this.keySet=fileCon.getClassificationBlock("Taxonomy").getKeySet();
			for(int key : keySet){
				map.put(key,fileCon.getClassificationBlock("Taxonomy").getSum(key));
			}
			this.assignmentMap = map;
		}catch(IOException io){
			io.printStackTrace();
		}
	}
}
