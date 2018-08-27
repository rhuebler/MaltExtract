package DatabaseAnalyzer;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Map;
import java.util.Set;
import java.util.logging.Level;
import java.util.logging.Logger;
import NCBI_MapReader.NCBI_MapReader;
import RMAAlignment.NOAOR;
import RMAAlignment.NOAORComparator;
import megan.rma6.ClassificationBlockRMA6;
import megan.rma6.RMA6File;

public class ReadDatabaseAnalyzer {
	private String inDir;
	private String fileName;
	private Set<Integer> keySet;
	private Set<Integer> allKeys;
	private Integer readCount;
	private Map<Integer,Integer> assignmentMap;
	private String output;
	private Logger log;
	private Logger warning;
	private NCBI_MapReader mapReader;

	// constructor
	public ReadDatabaseAnalyzer(String inDir, String name,
			Logger log, Logger warning, NCBI_MapReader reader){
		this.inDir = inDir;
		this.fileName =  name;
		this.log = log;
		this.warning = warning;
		this.mapReader = reader;
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
	public String getOutput() {
		return output;
	}
	protected String getName(int taxId){
		String name;
		if(mapReader.getNcbiIdToNameMap().get(taxId) != null)
			name = mapReader.getNcbiIdToNameMap().get(taxId).replace(' ', '_').replace('\'', '_').replace('#', '_').replace('$', '_');
		else if(taxId == 0)
			name="NA";
		else
			name = "unassignedName";
		return name;
	}
	// process
	private void process(){
		try{
			log.log(Level.INFO, "Scanning File: "+ fileName);
			ArrayList<NOAOR> noarList= new ArrayList<NOAOR>();
			RMA6File rma6File = new RMA6File(inDir+fileName, "r");
			Long location = rma6File.getFooterSectionRMA6().getStartClassification("Taxonomy");
		    if (location != null) {
		        ClassificationBlockRMA6 cl = new ClassificationBlockRMA6("Taxonomy");
		        cl.read(location, rma6File.getReader());
		      //if no taxa list provided use all taxa 
		        this.keySet = cl.getKeySet();
		        this.allKeys = cl.getKeySet();
		        this.readCount = (int) rma6File.getFooterSectionRMA6().getNumberOfReads();// read in all counts
		        for(int key : keySet){
		        	if(allKeys.contains(key))
		        		noarList.add(new NOAOR(cl.getSum(key),null,key));
		        	else
		       			noarList.add(new NOAOR(0, null, key));
		        }
		    }else{
		    	warning.log(Level.SEVERE,fileName+" has no taxonomy block");
		    }
			rma6File.close();
			NOAORComparator comparator = new NOAORComparator();
			noarList.sort(comparator);
			String output = fileName+"\t";
			for(int i =0; i< 10; i++) {
				if(i<noarList.size()) {
					output += getName(noarList.get(i).getTaxID())+";"+noarList.get(i).getSize()+"\t";
				}else{
					output += "NA;NA\t";
				}
			}
			this.output = output;
		}catch(IOException io){
			warning.log(Level.SEVERE,"Cannot open File",io);
		}
	}
}
