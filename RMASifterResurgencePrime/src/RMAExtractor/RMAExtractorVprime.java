package RMAExtractor;

//TODO benchmark on more files and threads on single thread version and on concurrent version to see whether or not output is consistent
import java.io.File;
import java.io.FilenameFilter;
import java.util.ArrayList;
import java.util.List;
import java.util.Scanner;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.ThreadPoolExecutor;

import NCBI_MapReader.NCBI_MapReader;
import NCBI_MapReader.NCBI_TreeReader;
import RMA6Processor.ConcurrentRMA6Processor;
import RMA6Processor.RMA6Processor;
import SummaryWriter.SummaryWriter;

public class RMAExtractorVprime {
	static FilenameFilter RMAFilter = new FilenameFilter() {
	    public boolean accept(File file, String name) {
	        if (name.endsWith(".rma6")) {
	            return true;
	        } else {
	            return false;
	        }
	    }
	};// ToDo make parallel
	//implement taxname.file reader 
	private static double topPercent = 0.01; // initialize with standard value;	
	private static List<String> fileNames = new ArrayList<String>();
	private static List<String> taxNames = new ArrayList<String>();
	private static String outDir;
	private static String inDir;	
	private static int numThreads = 1;
	private static ThreadPoolExecutor executor;
	private static void destroy(){
		executor.shutdown();
	}
	public static void main(String[] args) throws Exception {
		System.out.println("Processing Input Paramenters");	
		inDir = args[0];
		File dir = new File(inDir);
		if(dir.isDirectory())
			{
			for(String name :dir.list(RMAFilter)){
				fileNames.add(name);}
		}else if(dir.isFile()&& dir.canRead()){
			fileNames.add(dir.getName());
			inDir=dir.getAbsolutePath().substring(0, dir.getAbsolutePath().indexOf(dir.getName()));
			}
		else{
			System.err.println("No RMA6 input File detected");
			System.err.println("Expected Parameter order:");
			System.err.println("-path/to/RMA6_files or RMA6_file (-outdir) -MeganID1 -MeganID2 -... (-topPercent)");
		}
		if(args[1].substring(1).matches("[\\w\\d/]+")){
			outDir=args[1];}
		else{
			System.out.println("Use input directory as output directory!");
			outDir=inDir;
		}
		for( String arg : args ){
			if(arg.endsWith("taxons.txt")){// should now match species Names but not threads species Names at least as more than 100 threads seem unreasonable characters long
				Scanner in = new Scanner(new File(arg));
				while(in.hasNext()){
					taxNames.add(in.nextLine().trim());
				}
				in.close();
				}else if(arg.matches("0\\.\\d+"))
				{
				topPercent=Double.parseDouble(arg);
				System.out.println("Custom topPercent parameter detected: "+topPercent);
				}else if (arg.matches("\\d+")){
				numThreads = Integer.parseInt(arg);
			}
		}
		if(numThreads >Runtime.getRuntime().availableProcessors())
			numThreads = Runtime.getRuntime().availableProcessors(); // enforce that not more processors than available to jvm should be used 
		System.out.println("Using " + numThreads +" threads");
		NCBI_MapReader mapReader = new NCBI_MapReader();// shared read access
		NCBI_TreeReader treeReader = new NCBI_TreeReader();// shared read access
		new File(outDir).mkdirs(); // create output directory if directory does not already exist
	    // iterate over files
	    
    	List<Integer> taxIDs= new  ArrayList<Integer>();
    	List<Future<RMA6Processor>> processedFiles = new ArrayList<>();
    	for(String name : taxNames){
    		if(mapReader.getNcbiNameToIdMap().get(name) != null)// catch if there is a mistake
    			taxIDs.add(mapReader.getNcbiNameToIdMap().get(name));
    		else
    			System.err.println(name + " has no assigned taxID and cannot be processed!");
    	}
    	
    	executor=(ThreadPoolExecutor) Executors.newFixedThreadPool(numThreads);
	    for(String fileName : fileNames){
	    	ConcurrentRMA6Processor task = new ConcurrentRMA6Processor(inDir, fileName, outDir, mapReader, treeReader,taxIDs, topPercent); 
	    	Future<RMA6Processor> future=executor.submit(task);
	    	processedFiles.add(future);
	    }//fileNames;
	    // wait for all threads to finish here  synchronize system resources but how?
	    destroy();

	  SummaryWriter sumWriter = new SummaryWriter(processedFiles,mapReader,outDir); 
	  sumWriter.writeSummary();
	}//main
}//class
