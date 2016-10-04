package utility;


import java.io.File;
import java.io.FileNotFoundException;
import java.io.FilenameFilter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Scanner;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.regex.Pattern;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;

import behaviour.Filter;
import behaviour.Taxas;

/**
 * This class is used To Parse Input Parameters for RMA Extractor to make input more flexible and less error prone
 * uses the console parameters as input and fill up all Parameter slots to control subsequent functions
 * @author huebler
 *
 */
public class InputParameterProcessor {
	/**
	 * @param String[] args all comandline parameters
	 * @throws none thrown all caught
	 * @return return values and parameters to run and control the program
	 */
	private double topPercent = 0.01; // initialize with standard value;	
	private List<String> fileNames = new ArrayList<String>();
	private List<String> taxNames = new ArrayList<String>();
	private String outDir;
	private int numThreads = 1;
	private int maxLength = 0;
	private double minPIdent = 0;
	private Filter behave = Filter.SCAN;
	private Taxas taxas = Taxas.ALL;
	private String tree_Path = "/projects1/clusterhomes/huebler/RMASifter/RMA_Extractor_Resources/";
	private boolean verbose = false;
	private Logger log;
	private Logger warning;
	private boolean reads = false;
	// constructor
	public InputParameterProcessor(String[] params ,Logger log, Logger warning){
		this.log = log;
		this.warning =  warning;
		process(params);	
	}
	// getters
	public boolean getBlastHits(){
		return reads;
	}
	public double getMinPIdent(){
		return this.minPIdent;
	}
	public List<String> getTaxNames(){
		return this.taxNames;
	}
	
	public List<String> getFileNames(){
		return this.fileNames;
	}
	public double getTopPercent(){
		return this.topPercent;
	}
	public String getOutDir(){
		return this.outDir;
	}
	public Filter getFilter(){
		return this.behave;
	}
	public Taxas getTaxas(){
		return this.taxas;
	}
	public int getNumThreads(){
		return this.numThreads;
	}
	public int getMaxLength(){
		return this.maxLength;
	}
	public String getTreePath(){
		return this.tree_Path;
	}
	public boolean isVerbose(){
		return this.verbose;
	}
    private void process(String[] parameters)
    {	FilenameFilter RMAFilter = new FilenameFilter() {
	    	public boolean accept(File file, String name) {
	    		if (name.endsWith(".rma6")) {
	    			return true;
	    		} else {
	    			return false;
	    		}
	    	}
    	};
    	 CommandLine commandLine;
    	 	// Short Flags
    	 	// edit distance 
    	    Option option_Input = Option.builder("i").longOpt("input").argName("Path/to/inDir or RMA6Files").hasArgs().desc("Input Directory or files").build();
    	    Option option_Output = Option.builder("o").longOpt("output").argName("Path/to/outDir").hasArg().desc("Output Directory").build();
    	    Option option_Taxons = Option.builder("t").longOpt("taxons").argName("Path/to/taxFile or Taxon in \"\"").hasArg().desc("Taxons to look up").build();
    	    
    	    // long flags
    	    Option option_Threads = Option.builder("p").longOpt("threads").argName("1..maxNumberOfCores").hasArg().optionalArg(true).desc("Number of Cores to run on").build();		
    	    Option option_TopPercent = Option.builder("a").longOpt("top").argName("0.0-0.99").hasArg().optionalArg(true).desc("Top Percent of Matches to Consider").build();
    	    Option option_Filter = Option.builder("f").longOpt("filter").argName("non,ancient,nonduplicate, scan").optionalArg(true).hasArg().desc("Specify the behaviour for run eg ancient").build();
    	    Option option_MaxLength = Option.builder().longOpt("maxLength").argName("maxLength").hasArg().optionalArg(true).desc("Set Maximum Read Length").build();
    	    Option option_minPercentIdent = Option.builder().longOpt("minPI").argName("minPI").hasArg().optionalArg(true).desc("Set Minumum Percent Identity").build(); 
    	    Option option_Help = Option.builder("h").longOpt("help").optionalArg(true).desc("Print Usage and shutdown").build();
    	    Option option_Path = Option.builder("r").longOpt("resources").hasArg().optionalArg(true).desc("Path to NCBI tre and map File").build();
    	    Option option_Verbose = Option.builder("v").longOpt("verbose").optionalArg(true).desc("How much output should be printed to screen").build();
    	    Option option_Reads = Option.builder().longOpt("reads").optionalArg(true).desc("Output Blas Hits when filtering for ancient Reads").build();
    	    Options options = new Options();
    	    CommandLineParser parser = new DefaultParser();

    	    options.addOption(option_Input);
    	    options.addOption(option_Output);
    	    options.addOption(option_Taxons);
    	    
    	    options.addOption(option_Threads);
    	    options.addOption(option_TopPercent);
    	    options.addOption(option_Filter);
    	    options.addOption(option_MaxLength);
    	    options.addOption(option_minPercentIdent);
    	    
    	    options.addOption(option_Help);
    	    options.addOption(option_Path);
    	    options.addOption(option_Verbose);
    	    options.addOption(option_Reads);

    	    try
    	    {
    	        commandLine = parser.parse(options, parameters);

    	        if (commandLine.hasOption("input"))
    	        {
    	        	log.log(Level.INFO,"Input Set to: ");
    	            for(String arg :commandLine.getOptionValues("input")){
    	            	File inFile = new File(arg);
    	            	try{
    	            		if(inFile.getCanonicalFile().isDirectory()){
    	            			 log.info(arg);
    	            			for(String name : inFile.list(RMAFilter))
    	            				this.fileNames.add(arg + name);
    	            		}else if(inFile.isFile()){
    	            			log.info(arg);
    	            			this.fileNames.add(arg);
    	            		}
    	            	}catch(IOException io){
    	            		warning.log(Level.SEVERE,"Can't open File", io);
    	            	}	
    	            }
    	        }

    	        if (commandLine.hasOption("output"))
    	        {
    	        	log.log(Level.INFO,"Output Directory set to: "+commandLine.getOptionValue("output"));
    	            this.outDir = commandLine.getOptionValue("output");
    	        }
    	        
    	        if (commandLine.hasOption("taxons"))
    	        {	this.taxas = Taxas.USER;
    	            for(String tax : commandLine.getOptionValues("taxons")){
    	            	log.info("Taxons File: ");
    	            	log.info(tax);
    	        	   File f = new File(tax);
    	        	   if(f.isFile()){
    	        		   try {
    	        			   Scanner	in = new Scanner(f);
    	        			   while(in.hasNext()){
    	        				   taxNames.add(in.nextLine().trim());
    	        			   }
    	        			   in.close();
    	        		   }catch (FileNotFoundException e) {
    	        		   e.printStackTrace();
    	        		   }
    	        	   }else{
    	        		   log.info("Added Taxon: ");
    	        		   log.info(tax +" to analysis");
    	        		   taxNames.add(tax); 
    	        	   }
    	           }	   
    			}
    					        
    	        
    	        if (commandLine.hasOption("threads"))
    	        {	
    	        	log.info("Trying to use "+commandLine.getOptionValue("threads")+" threads");
    	            this.numThreads=Integer.parseInt(commandLine.getOptionValue("threads"));
    	            if(this.numThreads > Runtime.getRuntime().availableProcessors()){
    	    			this.numThreads = Runtime.getRuntime().availableProcessors(); // enforce that not more processors than available to jvm should be used 
    	    			log.warning("The Number of Threads higher than than Number available to System");
    	            }
    	    		log.log(Level.INFO,"Using " + numThreads +" threads");
    	        }

    	        if (commandLine.hasOption("filter"))
    	        {
    	            if(Pattern.compile(Pattern.quote("non"), Pattern.CASE_INSENSITIVE).matcher(commandLine.getOptionValue("filter")).find()){
    	            	this.behave = Filter.NON;
    	            	log.log(Level.INFO,"Custom Behaviour set to: "+commandLine.getOptionValue("filter"));
    	            }else if(Pattern.compile(Pattern.quote("ancient"), Pattern.CASE_INSENSITIVE).matcher(commandLine.getOptionValue("filter")).find()){
    	            	this.behave = Filter.ANCIENT;
    	            	log.log(Level.INFO,"Custom Behaviour set to: "+commandLine.getOptionValue("filter"));
    	            }else if(Pattern.compile(Pattern.quote("nonduplicate"), Pattern.CASE_INSENSITIVE).matcher(commandLine.getOptionValue("filter")).find()){
    	            	this.behave = Filter.NONDUPLICATES;
    	            	log.log(Level.INFO,"Custom Behaviour set to: "+commandLine.getOptionValue("filter"));
    	            }else if(Pattern.compile(Pattern.quote("scan"), Pattern.CASE_INSENSITIVE).matcher(commandLine.getOptionValue("filter")).find()){
    	            	this.behave = Filter.SCAN;
    	            	log.log(Level.INFO,"Custom Behaviour set to: "+commandLine.getOptionValue("filter"));
    	            }
    	            else if(Pattern.compile(Pattern.quote("all"), Pattern.CASE_INSENSITIVE).matcher(commandLine.getOptionValue("filter")).find()){
    	            	this.behave = Filter.ALL;
    	            	log.log(Level.INFO,"Custom Behaviour set to: "+commandLine.getOptionValue("filter"));
    	            }
    	            
    	        }

    	        if (commandLine.hasOption("top"))
    	        {
    	        	log.log(Level.INFO,"Top Percent value set to: "+commandLine.getOptionValue("top"));
    	            this.topPercent = Double.parseDouble(commandLine.getOptionValue("top"));
    	        }
    	        
    	        if (commandLine.hasOption("maxReadLength"))
    	        {
    	        	log.log(Level.INFO,"maximum Read Length set to: " + commandLine.getOptionValue("maxReadLength"));
    	            this.maxLength = Integer.parseInt(commandLine.getOptionValue("maxReadLength"));
    	        }
    	        
    	        if(commandLine.hasOption("minPI")){
    	        	this.minPIdent = Double.parseDouble(commandLine.getOptionValue("minPI"));
    	        }
    	        
    	        if(commandLine.hasOption("resources")){
    	        	File f = new File(commandLine.getOptionValue("resources"));
    	        	if(f.isDirectory()){
    	        		this.tree_Path = commandLine.getOptionValue("resources");
    	        	}
    	        }
    	        if(commandLine.hasOption('v')){
    	        	this.verbose = true;
    	        }
    	        if((commandLine.hasOption("reads") && behave == Filter.ANCIENT)||
    	        	(commandLine.hasOption("reads") && behave == Filter.ALL)	){
    	        	this.reads = true;
    	        }
    	        if(commandLine.hasOption("h")){
    	        	String header = "RMAExtractor concurrent alpha";
    	    	    String footer = "In case you encounter an error drop an email to huebler@shh.mpg.de with useful description";
    	    	    HelpFormatter formatter = new HelpFormatter();
    	    	    formatter.printHelp("RMAExtractor", header, options, footer, true);   
    	    	    System.exit(0);
    	        }
    	        if(!commandLine.hasOption("h")  && (!commandLine.hasOption("output") || !commandLine.hasOption("input"))){
    	        	warning.log(Level.SEVERE,"Please, specifiy input files or input directories and output directory");
    	        	System.exit(1);
    	        }
    	        String[] remainder = commandLine.getArgs();
    	        if(remainder.length != 0){
    	           warning.log(Level.SEVERE,"Remaining arguments: ");
    	            for (String argument : remainder){
    	              warning.log(Level.WARNING,argument);
    	            }
    	          }
    	       
    	    }
    	    catch (ParseException exception){
    	    	warning.log(Level.SEVERE, "Parse error: " + exception);
    	    }
    }
  
}