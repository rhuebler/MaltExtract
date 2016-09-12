package RMAExtractor;

//TODO benchmark on more files and threads on single thread version and on concurrent version to see whether or not output is consistent
//TODO adress errors by try catch if no other way 
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.ThreadPoolExecutor;

import NCBI_MapReader.NCBI_MapReader;
import NCBI_MapReader.NCBI_TreeReader;
import RMA6Processor.ConcurrentRMA6Processor;
import RMA6Processor.ConcurrentRMA6Scanner;
import RMA6Processor.RMA6Processor;
import RMA6Processor.RMA6Scanner;
import behaviour.Filter;
import behaviour.Taxas;
import utility.InputParameterProcessor;
import utility.ScanSummaryWriter;
import utility.SummaryWriter;
/**
 * Essentially the Main Class of RMA Extractor is a concurrent Program
 * At start up passes all command line arguments into InputProcessor
 * Input Processor grabs all arguments from the specified flags
 * and default values to unspecified arguments
 * than initializes instances of NCBI_Map_Reader and NCBI_TreeReader 
 * it creates Instances of RMA6Processor a class designed to Process one RMA6File 
 * and passes future instances of these to summary writer to write summary output at the end
 * @author huebler
 *
 */
public class RMAExtractor {
	
	private static ThreadPoolExecutor executor;
	
	private static void destroy(){
		executor.shutdown();
	}
	public static void main(String[] args) throws IOException  {
		long startTime = System.nanoTime();
		InputParameterProcessor inProcessor = new InputParameterProcessor(args);
		NCBI_MapReader mapReader = new NCBI_MapReader(inProcessor.getTreePath());
		new File(inProcessor.getOutDir()).mkdirs();
		new File(inProcessor.getOutDir()+"/readDist/").mkdirs();
		if(inProcessor.wantReadInf()){
			new File(inProcessor.getOutDir()+"/editDistance/").mkdirs();
			new File(inProcessor.getOutDir()+"/percentIdentity/").mkdirs();
		}
		List<Integer> taxIDs= new  ArrayList<Integer>();
		if(inProcessor.getTaxas() == Taxas.USER){
			for(String name : inProcessor.getTaxNames()){
				if(mapReader.getNcbiNameToIdMap().get(name) != null)// catch if there is a mistake
					taxIDs.add(mapReader.getNcbiNameToIdMap().get(name));
				else
					System.err.println(name + " has no assigned taxID and cannot be processed!");
			}
    	}
		
		executor=(ThreadPoolExecutor) Executors.newFixedThreadPool(inProcessor.getNumThreads());//intialize concurrent thread executor 
		
		if(inProcessor.getFilter() != Filter.SCAN){
			List<Future<RMA6Processor>> processedFiles = new ArrayList<>();
			 NCBI_TreeReader treeReader = new NCBI_TreeReader(inProcessor.getTreePath());
    		for(String fileName : inProcessor.getFileNames()){
    			File f = new File(fileName);
    			 // every tree has its own copy of this now to avoid concurrency issues
    			ConcurrentRMA6Processor task = new ConcurrentRMA6Processor(f.getParentFile().getCanonicalFile() + "/", f.getName(), inProcessor.getOutDir(), 
	    			mapReader, treeReader,taxIDs, inProcessor.getTopPercent(),inProcessor.getMaxLength(),inProcessor.getFilter(), 
	    			inProcessor.getTaxas(), inProcessor.wantReadInf(), inProcessor.isVerbose());
    			Future<RMA6Processor> future=executor.submit(task);
    			processedFiles.add(future);
    			
    		}//fileNames;
	    // wait for all threads to finish here currently no concurrency errors or deadlocks but this would be the place where it would fall apart 
	    destroy();
	    SummaryWriter sumWriter = new SummaryWriter(processedFiles,mapReader,inProcessor.getOutDir()); 
	    sumWriter.writeSummary();
	  }else{// TODO add functionality to support abstract file paths 
		  List<Future<RMA6Scanner>> scannerList = new ArrayList<Future<RMA6Scanner>>();
		  NCBI_TreeReader treeReader = new NCBI_TreeReader(inProcessor.getTreePath());
		  // every tree has its own copy of this now to avoid concurrency issues
		  for(String fileName : inProcessor.getFileNames()){
			 
			 File f = new File(fileName);
			 ConcurrentRMA6Scanner task = new ConcurrentRMA6Scanner(f.getParent()+"/", f.getName(),inProcessor.getTaxas(),taxIDs, treeReader);
			 Future<RMA6Scanner> future = executor.submit(task);
			 scannerList.add(future);
		  }
		  destroy();
		  ScanSummaryWriter writer = new ScanSummaryWriter(scannerList, mapReader);
		  writer.write(inProcessor.getOutDir());
	  }
		long endTime = System.nanoTime();
		System.out.println("Runtime: "+ (endTime - startTime)/1000000000 +" Seconds");
	}//main
	
}//class
