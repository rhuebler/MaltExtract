package utility;


import java.io.File;
import java.io.FileNotFoundException;
import java.io.FilenameFilter;
import java.util.ArrayList;
import java.util.List;
import java.util.Scanner;

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
	private double topPercent = 0.01; // initialize with standard value;	
	private List<String> fileNames = new ArrayList<String>();
	private List<String> taxNames = new ArrayList<String>();
	private String outDir;
	private int numThreads = 1;
	private int maxLength = 0;
	private Filter behave = Filter.NON;
	private Taxas taxas = Taxas.ALL;
	private String tree_Path = "/projects1/clusterhomes/huebler/RMASifter/RMA_Extractor_Resources/";
	// constructor
	public InputParameterProcessor(String[] params){
		process(params);
	}
	// getters
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
    	    Option option_Input = Option.builder("input").argName("Path/to/inDir or RMA6Files").hasArgs().required().desc("Input Directory or file").build();
    	    Option option_Output = Option.builder().longOpt("output").argName("Path/to/outDir").hasArg().required().desc("Output Directory").build();
    	    Option option_Taxons = Option.builder("taxons").argName("Path/to/taxFile or Taxon in \"\"").hasArg().desc("File with taxons to look up").build();
    	    
    	    // long flags
    	    Option option_Threads = Option.builder().longOpt("threads").argName("1..maxNumberOfCores").hasArg().optionalArg(true).desc("Number of Cores to run on").build();		
    	    Option option_TopPercent = Option.builder().longOpt("top").argName("0.0-0.99").hasArg().optionalArg(true).desc("Top Percent of Matches to Consider").build();
    	    Option option_Filter = Option.builder().longOpt("filter").argName("non,ancient,nonduplicate, scan").optionalArg(true).hasArg().desc("Specify the behaviour for run eg ancient").build();
    	    Option option_MaxLength = Option.builder().longOpt("maxReadLength").argName("maxLength").hasArg().optionalArg(true).desc("Set Maximum ReadLength").build();
    	    Option option_Help = Option.builder("h").optionalArg(true).desc("Print Usage and shutdown").build();
    	    Option option_Path = Option.builder().longOpt("tree").hasArg().optionalArg(true).desc("Path to NCBI tre and map File").build();
    	    Options options = new Options();
    	    CommandLineParser parser = new DefaultParser();

    	    options.addOption(option_Input);
    	    options.addOption(option_Output);
    	    options.addOption(option_Taxons);
    	    
    	    options.addOption(option_Threads);
    	    options.addOption(option_TopPercent);
    	    options.addOption(option_Filter);
    	    options.addOption(option_MaxLength);
    	    options.addOption(option_Help);
    	    options.addOption(option_Path);
    	    

    	    try
    	    {
    	        commandLine = parser.parse(options, parameters);

    	        if (commandLine.hasOption("input"))
    	        {
    	            System.out.print("Input Set to: ");
    	            for(String arg :commandLine.getOptionValues("input")){
    	            	File inFile = new File(arg);
    	            	if(inFile.isDirectory()){
    	            		  System.out.println(arg);
    	            		for(String name : inFile.list(RMAFilter))
    	            		this.fileNames.add(arg + name);
    	            	}else if(inFile.isFile()){
    	            		System.out.println(arg);
    	            		this.fileNames.add(arg);
    	            }
    	            }
    	        }

    	        if (commandLine.hasOption("output"))
    	        {
    	            System.out.print("Output Directory set to: ");
    	            System.out.println(commandLine.getOptionValue("output"));
    	            this.outDir = commandLine.getOptionValue("output");
    	        }
    	        
    	        if (commandLine.hasOption("taxons"))
    	        {	this.taxas = Taxas.USER;
    	            for(String tax : commandLine.getOptionValues("taxons")){
    	               System.out.println("Taxons File: ");
        	           System.out.println(tax);
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
    	        	   System.out.print("Added Taxon: ");
    	        	   System.out.println(tax +" to analysis");
    	        	   taxNames.add(tax); 
    	        	   }
    	           }	   
    			}
    					        
    	        
    	        if (commandLine.hasOption("threads"))
    	        {
    	            System.out.print("Trying to use ");
    	            System.out.println(commandLine.getOptionValue("threads")+" threads");
    	            this.numThreads=Integer.parseInt(commandLine.getOptionValue("threads"));
    	            if(this.numThreads > Runtime.getRuntime().availableProcessors()){
    	    			this.numThreads = Runtime.getRuntime().availableProcessors(); // enforce that not more processors than available to jvm should be used 
    	    			System.err.println("The Number of Threads higher than than Number available to System");
    	            }
    	    		System.out.println("Using " + numThreads +" threads");
    	        }

    	        if (commandLine.hasOption("filter"))
    	        {
    	            if(Pattern.compile(Pattern.quote("non"), Pattern.CASE_INSENSITIVE).matcher(commandLine.getOptionValue("filter")).find()){
    	            	this.behave = Filter.NON;
    	            	System.out.println("Custom Behaviour set to: ");
        	            System.out.println(commandLine.getOptionValue("filter"));
    	            }else if(Pattern.compile(Pattern.quote("ancient"), Pattern.CASE_INSENSITIVE).matcher(commandLine.getOptionValue("filter")).find()){
    	            	this.behave = Filter.ANCIENT;
    	            	System.out.println("Custom Behaviour set to: ");
        	            System.out.println(commandLine.getOptionValue("filter"));
    	            }else if(Pattern.compile(Pattern.quote("nonduplicate"), Pattern.CASE_INSENSITIVE).matcher(commandLine.getOptionValue("filter")).find()){
    	            	this.behave = Filter.NONDUPLICATES;
    	            	System.out.println("Custom Behaviour set to: ");
        	            System.out.println(commandLine.getOptionValue("filter"));
    	            }else if(Pattern.compile(Pattern.quote("scan"), Pattern.CASE_INSENSITIVE).matcher(commandLine.getOptionValue("filter")).find()){
    	            	this.behave = Filter.SCAN;
    	            	System.out.println("Custom Behaviour set to: ");
        	            System.out.println(commandLine.getOptionValue("filter"));
    	            }
    	            
    	        }

    	        if (commandLine.hasOption("top"))
    	        {
    	            System.out.println("Top Percent value set to: ");
    	            System.out.println(commandLine.getOptionValue("top"));
    	            this.topPercent = Double.parseDouble(commandLine.getOptionValue("top"));
    	        }
    	        
    	        if (commandLine.hasOption("maxReadLength"))
    	        {
    	            System.out.println("maximum Read Length set to: ");
    	            System.out.println(commandLine.getOptionValue("maxReadLength"));
    	            this.maxLength = Integer.parseInt(commandLine.getOptionValue("maxReadLength"));
    	        }
    	        
    	        if(commandLine.hasOption("tree")){
    	        	File f = new File(commandLine.getOptionValue("tree"));
    	        	if(f.isDirectory()){
    	        		this.tree_Path = commandLine.getOptionValue("tree");
    	        	}
    	        }
    	        
    	        if(commandLine.hasOption("h")){
    	        	String header = "RMAExtractor concurrent alpha";
    	    	    String footer = "In case you encounter an error drop an email to huebler@shh.mpg.de with useful description";
    	    	    HelpFormatter formatter = new HelpFormatter();
    	    	    formatter.printHelp("RMAExtractor", header, options, footer, true);   
    	    	    System.exit(0);
    	        }
    	        
    	        {
    	            String[] remainder = commandLine.getArgs();
    	            System.out.print("Remaining arguments: ");
    	            for (String argument : remainder)
    	            {
    	                System.out.print(argument);
    	                System.out.print(" ");
    	            }

    	            System.out.println();
    	        }
    	       
    	    }
    	    
    	    catch (ParseException exception)
    	    {
    	        System.out.print("Parse error: ");
    	        System.out.println(exception.getMessage());
    	    }
    }
  
}