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

import behaviour.Behaviour;

/**
 * This class is used To Parse Input Parameters for RMA Extractor to meke input more flexible and less error prone
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
	private int maxLength = 300;
	private Behaviour behave = Behaviour.ALL;
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
	public Behaviour getBehaviour(){
		return this.behave;
	}
	public int getNumThreads(){
		return this.numThreads;
	}
	public int getMaxLength(){
		return this.maxLength;
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
    	    Option option_Taxons = Option.builder("taxons").argName("Path/to/taxFile or Taxon in \"\"").hasArg().required().desc("File with taxons to look up").build();
    	    
    	    // long flags
    	    Option option_Threads = Option.builder().longOpt("threads").argName("1..maxNumberOfCores").hasArg().optionalArg(true).desc("Number of Cores to run on").build();		
    	    Option option_TopPercent = Option.builder().longOpt("top").argName("0.0-0.99").hasArg().desc("Top Percent of Matches to Consider").build();
    	    Option option_Behaviour = Option.builder().longOpt("behaviour").argName("all,ancient,nonduplicate, scan").optionalArg(true).hasArg().desc("Specify the behaviour for run eg ancient").build();
    	    Option option_MaxLength = Option.builder().longOpt("maxReadLength").argName("maxLength").optionalArg(true).hasArg().desc("Set Maximum ReadLength").build();
    	    Options options = new Options();
    	    CommandLineParser parser = new DefaultParser();

    	    options.addOption(option_Input);
    	    options.addOption(option_Output);
    	    options.addOption(option_Taxons);
    	    
    	    options.addOption(option_Threads);
    	    options.addOption(option_TopPercent);
    	    options.addOption(option_Behaviour);
    	    options.addOption(option_MaxLength);
    	    
    	    String header = "RMAExtractor concurrent alpha";
    	    String footer = "In case you encounter an error drop an email to huebler@shh.mpg.de with useful description";
    	    HelpFormatter formatter = new HelpFormatter();
    	    formatter.printHelp("RMAExtractor", header, options, footer, true);    

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
    	        {
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

    	        if (commandLine.hasOption("behaviour"))
    	        {
    	            if(Pattern.compile(Pattern.quote("all"), Pattern.CASE_INSENSITIVE).matcher(commandLine.getOptionValue("behaviour")).find()){
    	            	this.behave = Behaviour.ALL;
    	            	System.out.println("Custom Behaviour set to: ");
        	            System.out.println(commandLine.getOptionValue("behaviour"));
    	            }else if(Pattern.compile(Pattern.quote("ancient"), Pattern.CASE_INSENSITIVE).matcher(commandLine.getOptionValue("behaviour")).find()){
    	            	this.behave = Behaviour.ANCIENT;
    	            	System.out.println("Custom Behaviour set to: ");
        	            System.out.println(commandLine.getOptionValue("behaviour"));
    	            }else if(Pattern.compile(Pattern.quote("nonduplicate"), Pattern.CASE_INSENSITIVE).matcher(commandLine.getOptionValue("behaviour")).find()){
    	            	this.behave = Behaviour.NONDUPLICATES;
    	            	System.out.println("Custom Behaviour set to: ");
        	            System.out.println(commandLine.getOptionValue("behaviour"));
    	            }else if(Pattern.compile(Pattern.quote("scan"), Pattern.CASE_INSENSITIVE).matcher(commandLine.getOptionValue("behaviour")).find()){
    	            	this.behave = Behaviour.SCAN;
    	            	System.out.println("Custom Behaviour set to: ");
        	            System.out.println(commandLine.getOptionValue("behaviour"));
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