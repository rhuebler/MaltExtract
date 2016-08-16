package SummaryWriter;
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
import java.util.concurrent.ExecutionException;
import java.util.concurrent.Future;

import NCBI_MapReader.NCBI_MapReader;
import RMA6Processor.RMA6Processor;

public class SummaryWriter {
	private NCBI_MapReader mapReader;
	private List<Future<RMA6Processor>> processedFiles;
	private Set<Integer> processedIDs;
	private List<String> summary;
	private String outDir;
	
	public SummaryWriter(List<Future<RMA6Processor>> pFiles, NCBI_MapReader mReader, String oDir) throws IOException, InterruptedException, ExecutionException{
		this.processedFiles = pFiles;
		//this.processedIDs = pIDs;	
		this.mapReader = mReader;
		this.outDir = oDir;
		setProcessedIds();
		prepareOutput();
	}
	void setProcessedIds() throws InterruptedException, ExecutionException{
		Set<Integer> pIDs = new HashSet<>();
		for(Future<RMA6Processor> future : processedFiles)
			pIDs.addAll(future.get().getContainedIDs());
		this.processedIDs = pIDs;
	}
	private void prepareOutput() throws InterruptedException, ExecutionException{
		   List<String> summary = new ArrayList<String>();
		   String header ="Taxon"; // could and should be its own function 
		   boolean first = true;	   
		   for(Future<RMA6Processor> future : processedFiles){
			   RMA6Processor current = future.get();
			   header+="\t"+current.getfileName();
			   HashMap<Integer,Integer> fileResults = current.getSumLine();
			   if(first ==true){
			   for(int id : processedIDs){
				   String line = mapReader.getNcbiIdToNameMap().get(id);
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
				   }//for
			   }//else
		   }// for outer   
		   summary.add(0,header);
		   this.summary =  summary;
	}
	public void writeSummary() throws IOException{
		System.out.println("Writing Summary txt File");
		Path file = Paths.get(outDir+"overallSummary"+".txt");
		Files.write(file, summary, Charset.forName("UTF-8"));
		System.out.println("Summary Done!");
	}
}
