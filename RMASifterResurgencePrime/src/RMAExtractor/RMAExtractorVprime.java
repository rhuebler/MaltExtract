package RMAExtractor;
// TODO how do I make you multi threaded also test every thing works on more than one file 
// TODO how long is the shortest taxa Name 
// TODO great everything seems to be up and ready 
import java.io.File;
import java.io.FilenameFilter;
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

import NCBI_MapReader.NCBI_MapReader;
import NCBI_MapReader.NCBI_TreeReader;
import RMA6Processor.RMA6Processor;

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
	public static void writeSummary(List<String> summary) throws IOException{
		System.out.println("Writing Summary txt File");
		Path file = Paths.get(outDir+"overallSummary"+".txt");
		Files.write(file, summary, Charset.forName("UTF-8"));
		System.out.println("Summary Done!");
	}
	
	
	private static double topPercent = 0.01; // initialize with standard value;	
	private static List<String> fileNames = new ArrayList<String>();
	private static List<String> taxNames = new ArrayList<String>();
	private static String outDir;
	private static String inDir;	
	private static int threads = 1;
	
	public static void main(String[] args) throws IOException {
		System.out.println("Processing Input Paramenters");	
		inDir = args[0].substring(1);
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
		if(new File(args[1].substring(1)).isDirectory()){
			outDir=args[1].substring(1);}
		else{
			System.out.println("Use input directory as output directory!");
			outDir=inDir;
		}
		for( String arg : args ){
			if(arg.matches("-[\\w\\s\\d]{4,}")){// should now match species Names but not threads species Names at least as more than 100 threads seem unreasonable characters long
				taxNames.add(arg.substring(1));
				}else if(arg.matches("-0\\.\\d+"))
				{
				topPercent=Double.parseDouble(arg.substring(1));
				System.out.println("Custom topPercent parameter detected: "+topPercent);
				}else if (arg.matches("-\\d+")){
				threads = Integer.parseInt(arg.substring(1));
			}
		}
		if(threads >= 32)
			threads = 32; // only allow maximum of 32 threads /// how many cores do exist on cluster node?
		System.out.println("Using " +threads+" threads");
		NCBI_MapReader mapReader =new NCBI_MapReader();// shared read access
		NCBI_TreeReader treeReader = new NCBI_TreeReader();// shared read access
		new File(outDir).mkdirs(); // create output directory if directory does not already exist
	    // iterate over files
	    List<String> summary = new ArrayList<String>();
    	List<Integer>taxIDs= new  ArrayList<Integer>();
    	List<RMA6Processor>processedFiles = new  ArrayList<RMA6Processor>();
    	for(String name : taxNames)
    		taxIDs.add(mapReader.getNcbiNameToIdMap().get(name));
    	Set<Integer> processedIDs = new HashSet<Integer>();
	    for(String fileName : fileNames){// make multi threaded here
	    	RMA6Processor processor = new RMA6Processor(inDir, fileName, outDir, mapReader, treeReader);
	    	processor.process(taxIDs, topPercent);// loop through file
	    	processedIDs.addAll(processor.getContainedIDs());//TODO synchronize
	    	processedFiles.add(processor);//TODO synchronize
	    }//fileNames;
	    // wait for all threads to finish here  synchronize system resources but how?
	    
	   String header ="Taxon";
	   boolean first = true;
	   for(RMA6Processor current : processedFiles){
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
	   writeSummary(summary);
	    
	}//main
	
}//class
