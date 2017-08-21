package utility;
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

import behaviour.Filter;
/**
 * Write overall Summy file from processed RMA6Files
 * @author huebler
 *
 */
public class SummaryWriter {
	/**
	 * @param List<Future<RMA6Processor>> processedFiles, NCBI_MapReader mapReader, Set<Integer> processedIDs,
	 * String outDir Logger,Warning
	 * @throws none thrown all caught
	 */
	private NCBI_MapReader mapReader;
	private Set<Integer> processedIDs;
	private String outDir;
	private Logger warning;
	private Filter behave;
	private HashMap<String, HashMap<Integer,Integer>> overallSum;
	private HashMap<String, HashMap<Integer,Integer>> ancientSum;
	private HashMap<String,Integer> totalCounts;
	public SummaryWriter(HashMap<String, HashMap<Integer,Integer>> overallSum, HashMap<String, HashMap<Integer,Integer>> ancientSum, HashMap<String,Integer> totalCounts, 
			HashSet<Integer> containedIDs, NCBI_MapReader mReader, String oDir, Logger warning, Filter behave){
		this.overallSum = overallSum;
		this.ancientSum = ancientSum;
		this.totalCounts = totalCounts;
		this.mapReader = mReader;
		this.outDir = oDir;
		this.warning = warning;
		this.behave = behave;
		
	}
	public void process(){

		if(behave == Filter.NON_ANCIENT){
			prepareOutput(Filter.NON, overallSum);
			prepareOutput(Filter.NON_ANCIENT, ancientSum);
			prepareTotalReads(Filter.NON);
		}else if(behave == Filter.NON){
			prepareOutput(Filter.NON, overallSum);
			prepareTotalReads(Filter.NON);
		}else if(behave == Filter.NON_ANCIENT){
			prepareOutput(Filter.NON_ANCIENT, ancientSum);
			prepareTotalReads(Filter.NON_ANCIENT);
		}
	}
	private void prepareTotalReads(Filter switcher){
		ArrayList<String> totalReads = new  ArrayList<String>();
		String reads = "Total_Count";
		String header ="FileName";
		for(String fileName:totalCounts.keySet()){
			header += "\t"+fileName;
			reads += "\t" + totalCounts.get(fileName);
		}
		totalReads.add(header);
		totalReads.add(reads);
		 if(switcher==Filter.NON){
			 writeSummary(totalReads,outDir+"/default/"+"RunSummary"+".txt");
		 }else if(switcher==Filter.ANCIENT){
			 writeSummary(totalReads,outDir+"/ancient/"+"RunSummary"+".txt");
		 }
	}
	private void prepareOutput(Filter switcher, HashMap<String, HashMap<Integer,Integer>> sumlines) {
		   List<String> summary = new ArrayList<String>();
		   String header ="Node"; // could and should be its own function 
		  
		   boolean first = true;	   
		   for(String fileName:sumlines.keySet()){
				  header+="\t" + fileName;
				  HashMap<Integer,Integer> fileResults = sumlines.get(fileName);			  
				  if(first ==true){
					  for(int id : processedIDs){
						  String line;
						if( mapReader.getNcbiIdToNameMap().get(id) != null){
							line = mapReader.getNcbiIdToNameMap().get(id).replace(' ', '_');
						}else{
							line = "unasigned";
						}
						if(fileResults.containsKey(id)){
							line+= "\t"+fileResults.get(id);
							summary.add(line);
						}else{
							line+= "\t"+0;
							summary.add(line);
						}
					   first = false;
					   }//for
				   	}else{
					   int i = 0;
					   for(int id : processedIDs){ 
						   String line = summary.get(i);
						   if(fileResults.containsKey(id)){
							   line+= "\t"+fileResults.get(id);
							   summary.set(i,line);
						   }else{
							   line+= "\t"+0;
							   summary.set(i,line);
						   }
					  
						   i++;
					   }
				   	}
		   }
			   
		  
		   summary.sort(null);
		   summary.add(0,header);
		   if(switcher==Filter.NON){
			   writeSummary(summary,outDir+"/default/"+"RunSummary"+".txt");
		   }else if(switcher==Filter.ANCIENT){
			   writeSummary(summary,outDir+"/ancient/"+"RunSummary"+".txt");
		   }
	}

	private void writeSummary(List<String> summary,String outDir) {
		try{
		Path file = Paths.get(outDir);
		Files.write(file, summary, Charset.forName("UTF-8"));
		}catch(IOException io){
			warning.log(Level.SEVERE, "Error", io);
			System.exit(1);
		}
	}
	
}
