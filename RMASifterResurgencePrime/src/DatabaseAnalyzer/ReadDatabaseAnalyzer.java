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
	private NCBI_MapReader reader;

	// constructor
	public ReadDatabaseAnalyzer(String inDir, String name,
			Logger log, Logger warning, NCBI_MapReader reader){
		this.inDir = inDir;
		this.fileName =  name;
		this.log = log;
		this.warning = warning;
		this.reader = reader;
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
			String output = fileName;
			for(int i =0; i< 10; i++) {
				if(i<noarList.size()) {
					System.out.println(reader.getNcbiIdToNameMap().get(noarList.get(i).getTaxID())+"\t"+noarList.get(i).getSize());
					output += reader.getNcbiIdToNameMap().get(noarList.get(i).getTaxID())+";"+noarList.get(i).getSize()+"\t";
				}else{
					output += "NA;NA\t";
				}
			}
			System.out.println(output);
			this.output = output;
		}catch(IOException io){
			warning.log(Level.SEVERE,"Cannot open File",io);
		}
	}
}
