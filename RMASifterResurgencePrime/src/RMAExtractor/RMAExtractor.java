package RMAExtractor;

import java.io.File; //TODO rethink how to adress diffferent filters
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.ThreadPoolExecutor;
import java.util.logging.FileHandler;
import java.util.logging.Handler;
import java.util.logging.Level;
import java.util.logging.Logger;

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
	/**
	 * @param String[] args all comandline parameters
	 * @throws none thrown all caught
	 */
	private static final Logger log = Logger.getLogger(RMAExtractor.class.getName());
	private static final Logger warning = Logger.getLogger("Error");
	private static ThreadPoolExecutor executor;
	
	private static void destroy(){
		executor.shutdown();
	}
	public static void main(String[] args){
		long startTime = System.nanoTime();
		InputParameterProcessor inProcessor = new InputParameterProcessor(args ,log, warning);
		Handler handler = null;
		try {
			handler = new FileHandler(inProcessor.getOutDir()+"log.txt");
		} catch (SecurityException | IOException e) {
			e.printStackTrace();
		}
		log.addHandler(handler);
		
		Handler error = null;
		try {
			error = new FileHandler(inProcessor.getOutDir()+"error.txt");
		} catch (SecurityException | IOException e) {
			e.printStackTrace();
		}
		warning.addHandler(error);
		log.log(Level.INFO, "Setting up Taxon Name and Taxon ID maps");
		NCBI_MapReader mapReader = new NCBI_MapReader(inProcessor.getTreePath());
		new File(inProcessor.getOutDir()).mkdirs();
		new File(inProcessor.getOutDir()+"/readDist/").mkdirs(); //TODO could break potentially on Windows systems
		if(inProcessor.wantReadInf()){
			new File(inProcessor.getOutDir()+"/editDistance/").mkdirs();
			new File(inProcessor.getOutDir()+"/percentIdentity/").mkdirs();
		}
		List<Integer> taxIDs= new  ArrayList<Integer>();
		if(inProcessor.getTaxas() == Taxas.USER){
			for(String name : inProcessor.getTaxNames()){
				if(mapReader.getNcbiNameToIdMap().get(name) != null)// catch if there is a mistake
					taxIDs.add(mapReader.getNcbiNameToIdMap().get(name));
				else{
					warning.log(Level.SEVERE, name + " has no assigned taxID and cannot be processed!");
				}
			}
    	}
		
		executor=(ThreadPoolExecutor) Executors.newFixedThreadPool(inProcessor.getNumThreads());//intialize concurrent thread executor 
		log.log(Level.INFO, "Setting up Phylogenetic Tree");
		NCBI_TreeReader treeReader = new NCBI_TreeReader(inProcessor.getTreePath());
		if(inProcessor.getFilter() != Filter.SCAN){
			List<Future<RMA6Processor>> processedFiles = new ArrayList<>();
    		for(String fileName : inProcessor.getFileNames()){
    			try{
    				File f = new File(fileName);
    				ConcurrentRMA6Processor task = new ConcurrentRMA6Processor(f.getParentFile().getCanonicalFile() + "/", f.getName(), inProcessor.getOutDir(), 
    						mapReader, treeReader,taxIDs, inProcessor.getTopPercent(),inProcessor.getMaxLength(),inProcessor.getFilter(), 
    						inProcessor.getTaxas(), inProcessor.wantReadInf(), inProcessor.isVerbose(), log, warning);
    				Future<RMA6Processor> future=executor.submit(task);
    				processedFiles.add(future);
    				System.gc();
    			}catch(IOException io){
    				warning.log(Level.SEVERE,"File not found",io);
       				}
    		}//fileNames;
	    // wait for all threads to finish here currently no concurrency errors or deadlocks but this would be the place where it would fall apart 
	    destroy();
	    SummaryWriter sumWriter = new SummaryWriter(processedFiles,mapReader,inProcessor.getOutDir(), warning); 
	    log.log(Level.INFO, "Writing Summary File");
	    sumWriter.writeSummary();
	  }else{
		  List<Future<RMA6Scanner>> scannerList = new ArrayList<Future<RMA6Scanner>>();
		  // every tree has its own copy of this now to avoid concurrency issues
		  for(String fileName : inProcessor.getFileNames()){
			 try{
			 File f = new File(fileName);
			 ConcurrentRMA6Scanner task = new ConcurrentRMA6Scanner(f.getParentFile().getCanonicalFile()+"/",
					 f.getName(),inProcessor.getTaxas(),taxIDs, treeReader, log, warning);
			 Future<RMA6Scanner> future = executor.submit(task);
			 scannerList.add(future);
			 }catch(IOException io){
				 warning.log(Level.SEVERE,"File not found",io);
			 }
		  }
		  destroy();
		  ScanSummaryWriter writer = new ScanSummaryWriter(scannerList, mapReader, warning);
		  log.log(Level.INFO, "Writing Scan Summary File");
		  writer.write(inProcessor.getOutDir());
	  }
		long endTime = System.nanoTime();
		log.log(Level.INFO,"Runtime: "+ (endTime - startTime)/1000000000 +" Seconds");
	}//main
	
}//class
