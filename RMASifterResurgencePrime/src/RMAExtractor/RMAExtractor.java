package RMAExtractor;

//TODO benchmark on more files and threads on single thread version and on concurrent version to see whether or not output is consistent
//TODO adress errors by try catch if no other way 
import java.io.File;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.ThreadPoolExecutor;

import NCBI_MapReader.NCBI_MapReader;
import NCBI_MapReader.NCBI_TreeReader;
import RMA6Processor.ConcurrentRMA6Processor;
import RMA6Processor.RMA6Processor;
import utility.InputParameterProcessor;
import utility.SummaryWriter;
/**
 * Essentially the Main Class of RMA Extractor is a concurrent Programm
 * At start up passes all comandline arguments into InputProcessor
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
	public static void main(String[] args) throws Exception {
		InputParameterProcessor inProcessor = new InputParameterProcessor(args);
		NCBI_MapReader mapReader = new NCBI_MapReader();// shared read access
		NCBI_TreeReader treeReader = new NCBI_TreeReader();// shared read access
		new File(inProcessor.getOutDir()).mkdirs(); // create output directory if directory does not already exist
	    // iterate over files
	    
    	List<Integer> taxIDs= new  ArrayList<Integer>();
    	List<Future<RMA6Processor>> processedFiles = new ArrayList<>();
    	for(String name : inProcessor.getTaxNames()){
    		if(mapReader.getNcbiNameToIdMap().get(name) != null)// catch if there is a mistake
    			taxIDs.add(mapReader.getNcbiNameToIdMap().get(name));
    		else
    			System.err.println(name + " has no assigned taxID and cannot be processed!");
    	}
    	
    	executor=(ThreadPoolExecutor) Executors.newFixedThreadPool(inProcessor.getNumThreads());
	    for(String fileName : inProcessor.getFileNames()){
	    	ConcurrentRMA6Processor task = new ConcurrentRMA6Processor(inProcessor.getInDir(), fileName, inProcessor.getOutDir(), 
	    			mapReader, treeReader,taxIDs, inProcessor.getTopPercent(),inProcessor.getMaxLength(),inProcessor.getBehaviour()); 
	    	Future<RMA6Processor> future=executor.submit(task);
	    	processedFiles.add(future);
	    }//fileNames;
	    // wait for all threads to finish here currently no conuccrency errors or deadlocks 
	    destroy();

	  SummaryWriter sumWriter = new SummaryWriter(processedFiles,mapReader,inProcessor.getOutDir()); 
	  sumWriter.writeSummary();
	}//main
}//class
