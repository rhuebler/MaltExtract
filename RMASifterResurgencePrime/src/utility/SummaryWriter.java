package utility;
import java.io.IOException;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.Future;
import java.util.logging.Level;
import java.util.logging.Logger;



import NCBI_MapReader.NCBI_MapReader;
import RMA6Processor.RMA6Processor;
import behaviour.Filter;
/**
 * Write overall Summary file from processed RMA6Files for MaltExtract run.
 * Outputs into a table for every node in every file how many reads are still assigned after filtering
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
	private List<Future<RMA6Processor>> processedFiles;
	private Set<Integer> processedIDs;
	private String outDir;
	private Logger warning;
	private Filter behave;
	public SummaryWriter(List<Future<RMA6Processor>> pFiles, NCBI_MapReader mReader, String oDir, Logger warning, Filter behave){
		this.processedFiles = pFiles;
		this.mapReader = mReader;
		this.outDir = oDir;
		this.warning = warning;
		this.behave = behave;
		
	}
	public void process(){
		setProcessedIds();
		switch(behave) {
		default:
			prepareOutput(behave);
			break;
		case NON_ANCIENT:
			prepareOutput(Filter.NON);
			prepareOutput(Filter.ANCIENT);
			break;
		case SRNA:
			prepareOutput(Filter.NON);
			prepareOutput(Filter.ANCIENT);
			break;
		}
	}
	private void setProcessedIds(){
		HashSet<Integer> pIDs = new HashSet<Integer>();
		String fileName="";
		for(Future<RMA6Processor> future : processedFiles) {
			try{
			if(future.get().getContainedIDs()!=null) {
				pIDs.addAll(future.get().getContainedIDs());
				fileName = future.get().getfileName();
				}
			}catch(InterruptedException ie){
			warning.log(Level.SEVERE, "Error in File "+fileName, ie);
			}catch(ExecutionException ee){
			warning.log(Level.SEVERE, "Error in File "+fileName, ee);
			}
		}
		this.processedIDs = pIDs;
	}
	private void prepareOutput(Filter switcher) {
		 HashMap<Integer,String> summary = new  HashMap<Integer,String>();
		   String header ="Node"; // could and should be its own function 
		   String reads = "Total_Count";
		   boolean first = true;	   
		   for(Future<RMA6Processor> future : processedFiles){
			   RMA6Processor current;
			   try {
				   current = future.get();
				   header+="\t" + current.getfileName();
				   reads += "\t" + current.getTotalCount();
				   HashMap<Integer,Integer> fileResults = new  HashMap<Integer,Integer>();
				   if(switcher==Filter.NON)
					   fileResults = current.getSumLine();
				   else if(switcher==Filter.ANCIENT)
					   fileResults = current.getAncientLine();
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
							  
						   }else{
							   line+= "\t"+0;
						   }
						   summary.put(id,line);
					   first = false;
					   }//for
				   	}else{
					   for(int id : processedIDs){ 
						   String line = summary.get(id);
						   if(fileResults.containsKey(id)){
							   line+= "\t"+fileResults.get(id);
							  
						   }else{
							   line+= "\t"+0;
						   }
						   summary.replace(id, line);
					   }
				   	}
				  
			}catch (InterruptedException e) {
				warning.log(Level.SEVERE, "Interuption Error" , e);
			} catch (ExecutionException e) {
				warning.log(Level.SEVERE, "Execution Error", e);
			}
			   
		   }//for
		   ArrayList<String> output = new ArrayList<String>();
		   for(int id: summary.keySet())
				output.add(summary.get(id));
		   output.sort(null);
		   output.add(0,header);
		  
		   ArrayList<String> counts = new ArrayList<String>();
		   counts.add(header);
		   counts.add(reads);
		   
		   if(switcher==Filter.NON){
			   writeSummary(output,outDir+"/default/"+"RunSummary"+".txt");
			   writeSummary(counts,outDir+"/default/"+"TotalCount"+".txt");
		   }else if(switcher==Filter.ANCIENT){
			   writeSummary(output,outDir+"/ancient/"+"RunSummary"+".txt");
			   writeSummary(counts,outDir+"/ancient/"+"TotalCount"+".txt");
		   }
		   
	}

	private void writeSummary(List<String> summary,String outDir) {
		try{
		Path file = Paths.get(outDir);
		Files.write(file, summary, StandardCharsets.UTF_8);
		}catch(IOException io){
			warning.log(Level.SEVERE, "Error", io);
		}
	}
}
