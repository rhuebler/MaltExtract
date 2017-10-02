package RMAExtractor;
//TODO load tree and Map file direkte benutzbarkeit durch unkundigen users 
//TODO add ability to retireve output in SAM
//TODO add larger intervals
//TODO add %higher than X covered

import java.io.File; 
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.ThreadPoolExecutor;
import java.util.logging.FileHandler;
import java.util.logging.Handler;
import java.util.logging.Level;
import java.util.logging.Logger;

import NCBI_MapReader.NCBI_MapReader;
import NCBI_MapReader.NCBI_TreeReader;
import RMA6Processor.ConcurrentRMA6Crawler;
import RMA6Processor.ConcurrentRMA6Scanner;
import RMA6Processor.RMA6BlastCrawler;
import RMA6Processor.RMA6Scanner;
import RMA6Processor.parallelFileNodeProcessor;
import behaviour.Filter;
import behaviour.Taxas;
import jloda.util.PeakMemoryUsageMonitor;
import utility.DirectoryCreator;
import utility.InputParameterProcessor;
import utility.ScanSummaryWriter;
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
	 * @param String[] args all command line parameters
	 * @throws none thrown all caught
	 */
	private static final Logger log = Logger.getLogger(RMAExtractor.class.getName());
	private static final Logger warning = Logger.getLogger("Error");
	private static ThreadPoolExecutor executor;
	// check how to produce Megan File from RMAExtractor 
	private static void destroy(){
		executor.shutdown();
	}
	public static void main(String[] args){
		long startTime = System.nanoTime();
		 PeakMemoryUsageMonitor.start();
		InputParameterProcessor inProcessor = new InputParameterProcessor(args ,log, warning);
		new File(inProcessor.getOutDir()).mkdirs();// make outdir before log handlers
		Handler handler = null;
		//Initialize Output Handler and Error Handler
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
			 warning.log(Level.SEVERE,"Interuption",e);
		}
		warning.addHandler(error);
		
		//Set taxon and ID maps
		log.log(Level.INFO, "Setting up Taxon Name and Taxon ID maps");
		NCBI_MapReader mapReader = new NCBI_MapReader(inProcessor.getTreePath());
		DirectoryCreator dCreator = new DirectoryCreator();
		dCreator.process(inProcessor.getFilter(),inProcessor.getOutDir(),inProcessor.getBlastHits(),inProcessor.wantReads(), inProcessor.wantMeganSummaries());
		List<Integer> taxIDs= new  ArrayList<Integer>();
		
		NCBI_TreeReader treeReader = new NCBI_TreeReader(inProcessor.getTreePath());
		
		//read in taxa list and convert to numbers
		if(inProcessor.getTaxas() == Taxas.USER){
			for(String name : inProcessor.getTaxNames()){
				if(mapReader.getNcbiNameToIdMap().get(name) != null){// catch if there is a mistake
					taxIDs.add(mapReader.getNcbiNameToIdMap().get(name));
					if(inProcessor.isVerbose())
						log.log(Level.INFO, mapReader.getNcbiNameToIdMap().get(name)+" added to analysis");
				}else{
					warning.log(Level.SEVERE, name + " has no assigned taxID and cannot be processed!");
				}
			}
    	}
		
		//intialize  thread pool executor
		
		log.log(Level.INFO, "Setting up Phylogenetic Tree");
		
		// run normal mode if neither crawl nor scan are used
		if(inProcessor.getFilter() != Filter.SCAN  && inProcessor.getFilter() != Filter.CRAWL ){
    		parallelFileNodeProcessor fileProcessor = new parallelFileNodeProcessor(inProcessor,(ArrayList<Integer>) taxIDs,mapReader,treeReader, log, warning);
    		fileProcessor.process();
	    log.log(Level.INFO, "Writing Summary File");
	  }else{ 
			  if(inProcessor.getFilter() == Filter.SCAN && inProcessor.getFilter() != Filter.CRAWL){// run scan if crawl is not set
			  executor=(ThreadPoolExecutor) Executors.newFixedThreadPool(inProcessor.getNumThreads());//intialize concurrent thread executor 
			  List<Future<RMA6Scanner>> scannerList = new ArrayList<Future<RMA6Scanner>>();
			  // every tree has its own copy of this now to avoid concurrency issues
			  for(String fileName : inProcessor.getFileNames()){
				 File f = new File(fileName);
				 ConcurrentRMA6Scanner task = new ConcurrentRMA6Scanner(f.getParent()+"/",
						 f.getName(),inProcessor.getTaxas(),taxIDs, treeReader, log, warning,inProcessor.getOutDir(),inProcessor.wantMeganSummaries());
				 Future<RMA6Scanner> future = executor.submit(task);
				 scannerList.add(future);
			  }
			  destroy();
			  ScanSummaryWriter writer = new ScanSummaryWriter(scannerList, mapReader, warning);
			  log.log(Level.INFO, "Writing Scan Summary File");
			  writer.write(inProcessor.getOutDir());
	  }else if(inProcessor.getFilter() == Filter.CRAWL ){// run crawl filter
		  executor=(ThreadPoolExecutor) Executors.newFixedThreadPool(inProcessor.getNumThreads());//intialize concurrent thread executor 
		  for(String fileName : inProcessor.getFileNames()){
			  File f = new File(fileName);
			  log.log(Level.INFO, "Crawl for file " + fileName);
				 for(int taxID:taxIDs){
					 ConcurrentRMA6Crawler crawler = new ConcurrentRMA6Crawler(f.getParent()+"/",f.getName(),
							  mapReader.getNcbiIdToNameMap().get(taxID),
							  inProcessor.getOutDir(),mapReader, warning, treeReader, inProcessor.getFilter(),inProcessor.wantReads());
					  Future<RMA6BlastCrawler> future =  executor.submit(crawler);
					  try {
						future.get();
					} catch (InterruptedException e) {
						e.printStackTrace();
						System.exit(1);
					} catch (ExecutionException e) {
						e.printStackTrace();
						System.exit(1);
					}
				}
		  }	
	  } 
		  log.log(Level.INFO, "Crawling Done");
		  destroy();
	  }// get runtime 
		long endTime = System.nanoTime();
		log.log(Level.INFO,"Runtime: "+ (endTime - startTime)/1000000000 +" Seconds");
		log.log(Level.INFO,"Peak memory: " + PeakMemoryUsageMonitor.getPeakUsageString());
		log.log(Level.INFO,"Shutdown");
		System.exit(0);// shutdown
	}//main
	
}//class