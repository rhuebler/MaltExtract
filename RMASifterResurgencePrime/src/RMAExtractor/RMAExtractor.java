package RMAExtractor;
//TODO load tree and Map file direkte benutzbarkeit durch unkundigen users 
//TODO add ability to retireve output in SAM
//TODO add larger intervals
//TODO add %higher than X covered

import java.io.File; 
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
import java.util.logging.SimpleFormatter;

import NCBI_MapReader.NCBI_MapReader;
import NCBI_MapReader.NCBI_TreeReader;
import RMA6Processor.ConcurrentRMA6Processor;
import RMA6Processor.ConcurrentRMA6Scanner;
import RMA6Processor.RMA6BlastCrawler;
import RMA6Processor.RMA6Processor;
import RMA6Processor.RMA6Scanner;
import behaviour.Filter;
import behaviour.Taxas;
import jloda.util.PeakMemoryUsageMonitor;
import utility.DirectoryCreator;
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
	 * @param String[] args all command line parameters
	 * @throws none thrown all caught
	 */
	private static final Logger log = Logger.getLogger(RMAExtractor.class.getName());
	private static final Logger warning = Logger.getLogger(RMAExtractor.class.getName());
	private static ThreadPoolExecutor executor;
	private static void destroy(){
		executor.shutdown();
	}
	public static void main(String[] args){
		long startTime = System.nanoTime();
		 PeakMemoryUsageMonitor.start();
		InputParameterProcessor inProcessor = new InputParameterProcessor(args ,log, warning);
		new File(inProcessor.getOutDir()).mkdirs();// make outdir before log handlers
		Handler handler = null;
		SimpleFormatter sf = new SimpleFormatter();
		
		
		//Initialize Output Handler and Error Handler
		try {
			handler = new FileHandler(inProcessor.getOutDir()+"log.txt");
			handler.setFormatter(sf);
		} catch (SecurityException | IOException e) {
			e.printStackTrace();
		}
		log.addHandler(handler);
		Handler error = null;
		try {
			error = new FileHandler(inProcessor.getOutDir()+"error.txt");
			error.setFormatter(sf);
		} catch (SecurityException | IOException e) {
			 warning.log(Level.SEVERE,"Interuption",e);
		}
		warning.addHandler(error);
		//Set taxon and ID maps
		log.log(Level.INFO,inProcessor.GetParameters());
		
		log.log(Level.INFO, "Setting up Taxon Name and Taxon ID maps");
		NCBI_MapReader mapReader = new NCBI_MapReader(inProcessor.getTreePath());
		DirectoryCreator dCreator = new DirectoryCreator();
		dCreator.process(inProcessor.getFilter(),inProcessor.getOutDir(),inProcessor.getBlastHits(),inProcessor.wantReads(), inProcessor.wantMeganSummaries());
		ArrayList<Integer> taxIDs= new  ArrayList<Integer>();
		log.log(Level.INFO, "Setting up Phylogenetic Tree");
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
		
		
		
		// run normal mode if neither crawl nor scan are used

		if(inProcessor.getFilter() != Filter.SCAN  && inProcessor.getFilter() != Filter.CRAWL ){
			executor=(ThreadPoolExecutor) Executors.newFixedThreadPool(inProcessor.getNumThreads());//intialize concurrent thread executor 
			log.log(Level.INFO, "Using "+executor.getCorePoolSize()+"cores");
			List<Future<RMA6Processor>> processedFiles = new ArrayList<>();
    		for(String fileName : inProcessor.getFileNames()){
    			File f = new File(fileName);
					ConcurrentRMA6Processor task = new ConcurrentRMA6Processor(inProcessor,f.getParent() + "/",f.getName(), taxIDs, mapReader, treeReader, log, warning);
						Future<RMA6Processor> future=executor.submit(task);
						processedFiles.add(future);
    		}//fileNames;
	    // wait for all threads to finish here currently no concurrency errors or deadlocks but this would be the place where it would fall apart 
    	destroy();
	    SummaryWriter sumWriter = new SummaryWriter(processedFiles,mapReader,inProcessor.getOutDir(), warning,inProcessor.getFilter()); 
	    sumWriter.process();
	    log.log(Level.INFO, "Writing Summary File");
	  }else{ 
			  if(inProcessor.getFilter() == Filter.SCAN && inProcessor.getFilter() != Filter.CRAWL){// run scan if crawl is not set
			  executor=(ThreadPoolExecutor) Executors.newFixedThreadPool(inProcessor.getNumThreads());//intialize concurrent thread executor 
			  log.log(Level.INFO, "Using ",executor.getCorePoolSize());
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
		  for(String fileName : inProcessor.getFileNames()){
			  File f = new File(fileName);
			  log.log(Level.INFO, "Crawl for file " + fileName);
				 for(int taxID:taxIDs){
					 RMA6BlastCrawler crawler = new RMA6BlastCrawler(f.getParent()+"/",f.getName(),
							  mapReader.getNcbiIdToNameMap().get(taxID),
							  inProcessor,mapReader,log, warning, treeReader
							  );
					 crawler.process();
				}
		  }	
	  } 
		  log.log(Level.INFO, "Crawling Done");
	  }// get runtime 
		long endTime = System.nanoTime();
		log.log(Level.INFO,"Runtime: "+ (endTime - startTime)/1000000000 +" Seconds");
		log.log(Level.INFO,"Peak memory: " + PeakMemoryUsageMonitor.getPeakUsageString());
		log.log(Level.INFO,"Shutdown");
		System.exit(0);// shutdown
	}//main
	
}//class