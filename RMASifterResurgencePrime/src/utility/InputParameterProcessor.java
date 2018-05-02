package utility;


import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.nio.file.Files;
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
	// initialize class with standard values;	
	private double topPercent = 0.01; 
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
	private boolean alignment = false;
	private boolean reads = false;
	private double minComplexity = 0;
	private boolean wantMeganSummaries = false;
	private boolean destackOff = false;
	private boolean deDupOff=false;
	private boolean downsampling = true;
	// constructor

	public InputParameterProcessor(String[] params ,Logger log, Logger warning){
		 this.log = log;
		 this.warning =  warning;
		process(params);	
	}
	public String GetParameters() {
		String input ="";
		for(String name : fileNames)
			input+=name+"\b";
		String tax ="";
		for(String name : taxNames)
			tax+=name+"\b";
			String line="--input "+input+"\n"
				+"--taxa "+tax+"\n"
				+"--output "+outDir+"\n"
				+"--filter "+behave+"\n"
				+"--top "+topPercent+"\n"
				+"--maxLength "+maxLength+"\n"
				+"--minPI "+minPIdent+"\n"
				+"--resources "+tree_Path+"\n"
				+"--threads "+numThreads+"\n";
				if(verbose)
					line+="--verbose\n";
				if(alignment)
					line+="--matches\n";
				if(reads)
					line+="--reads\n";
				if(maxLength>0)
					line+="--maxReadLength "+maxLength+"\n";
				if(minComplexity>0)
					line+="--minComp "+minComplexity+"\n";
				if(wantMeganSummaries)
					line+="--wantMeganSummaries "+"\n";
				if(destackOff)
					line+="--destackingOff "+"\n";
				if(deDupOff)
					line+="--dupRemOff "+"\n";
				if(!downsampling)
					line+="--downSampOff "+"\n";
				return line;
	}
	// getters for parameters
	public boolean downsampling(){
		return downsampling;
	}
	public boolean turnDestackingOff(){
		return destackOff;
	}
	public double getMinComplexity(){
		return minComplexity;
	}
	public boolean wantMeganSummaries(){
		return  wantMeganSummaries;
	}
	public boolean getBlastHits(){
		return  alignment;
	}
	public double getMinPIdent(){
		return  minPIdent;
	}
	public List<String> getTaxNames(){
		return  taxNames;
	}
	
	public List<String> getFileNames(){
		return  fileNames;
	}
	public double getTopPercent(){
		return  topPercent;
	}
	public String getOutDir(){
		return  outDir;
	}
	public Filter getFilter(){
		return  behave;
	}
	public Taxas getTaxas(){
		return  taxas;
	}
	public int getNumThreads(){
		return  numThreads;
	}
	public int getMaxLength(){
		return  maxLength;
	}
	public String getTreePath(){
		return  tree_Path;
	}
	public boolean isVerbose(){
		return  verbose;
	}

	public boolean wantReads(){
		return  reads;
	}
	public boolean getDeDupOff(){
		return  deDupOff;
	}
	private void readTaxList(File f) throws IOException{
		try {
			 Scanner	in = new Scanner(f.getCanonicalFile());
			 while(in.hasNext()){
				 taxNames.add(in.nextLine().trim().replace('_', ' '));
			 }
			 in.close();
		 	}catch (FileNotFoundException e) {
		 warning.log(Level.WARNING,"File Not Found",e);
		}		
	}	
	private void readFileList(File f) throws IOException{
		try {
			 Scanner	in = new Scanner(f.getCanonicalFile());
			 while(in.hasNext()){
				 fileNames.add(in.nextLine());
			 }
			 in.close();
		 	}catch (FileNotFoundException e) {
		 warning.log(Level.WARNING,"File Not Found",e);
		}		
	}	
	private void process(String[] parameters){	
    	 CommandLine commandLine;
    	 	// Short Flags Are necessary parameters that are necessary for any run
    	    Option option_Input = Option.builder("i").longOpt("input").argName("Input").hasArgs().desc("Specify input directory or RMA6 files").build();
    	    Option option_Output = Option.builder("o").longOpt("output").argName("Outdir").hasArg().desc("Specify out directory").build();
    	    Option option_Taxons = Option.builder("t").longOpt("taxa").argName("Taxons").hasArgs().desc("Target species or List of targets").build();
    	    
    	    // long flags
    	    Option option_Threads = Option.builder("p").longOpt("threads").argName("Integer").hasArg().optionalArg(true).desc("How many cores to use?").build();		
    	    Option option_TopPercent = Option.builder("a").longOpt("top").argName("Double").hasArg().optionalArg(true).desc("Use top scoring 0.XX of alignments by defualt 0.01").build();
    	    Option option_Filter = Option.builder("f").longOpt("filter").argName("String").optionalArg(true).hasArg().desc("Use chosen filter full (def_anc), ancient, default, crawl, scan, srna").build();
    	    Option option_MaxLength = Option.builder().longOpt("maxReadLength").argName("Integer").hasArg().optionalArg(true).desc("Set maximum read length").build();
    	    Option option_minPercentIdent = Option.builder().longOpt("minPI").argName("Double").hasArg().optionalArg(true).desc("Set minimum percent identity to XX.X").build(); 
    	    Option option_Help = Option.builder("h").longOpt("help").optionalArg(true).desc("Print Help").build();
    	    Option option_Path = Option.builder("r").longOpt("resources").hasArg().argName("String").optionalArg(true).desc("Set path to required ncbi files").build();
    	    Option option_Verbose = Option.builder("v").longOpt("verbose").optionalArg(true).desc("How much output to print to screen").build();
    	    Option option_Alignment = Option.builder().longOpt("matches").optionalArg(true).desc("Retrieve Alignments").build();
    	    Option option_Reads = Option.builder().longOpt("reads").optionalArg(true).desc("Retrieve Reads").build();
    	   // Option option_Crawl = Option.builder().longOpt("crawl").optionalArg(true).desc("Use all alignments for damage and edit distance").build();
    	    Option option_minComplexity = Option.builder().longOpt("minComp").hasArg().argName("Double").optionalArg(true).desc("Use minimum complexity").build();
    	   // Option option_List = Option.builder().longOpt("list").hasArg().argName("list").optionalArg(true).desc("Decide on which build in list to use (not enabled yet)").build();
    	    Option option_MeganSummaries = Option.builder().longOpt("meganSummary").optionalArg(true).desc("Return Megan Summary Files").build();
    	    Option option_DeStackOff = Option.builder().longOpt("destackingOff").optionalArg(true).desc("Turn Off automated stacked Read Removal only useful in >1 coverage data").build();
    	    Option option_DeDupOff = Option.builder().longOpt("dupRemOff").optionalArg(true).desc("Turn Off automated pcr duplicate removal useful in >1 coverage data").build();
    	    Option option_Downsampling = Option.builder().longOpt("downSampOff").optionalArg(true).desc("Turn Off automatic downsampling on nodes with more than 10.000 assigned reads").build();
    	    Options options = new Options();
    	    
    	    // add all parameters to the parser
    	    CommandLineParser parser = new DefaultParser();

    	    options.addOption(option_Input);
    	    options.addOption(option_Output);
    	    options.addOption(option_Taxons);
    	    
    	    options.addOption(option_Threads);
    	    options.addOption(option_TopPercent);
    	    options.addOption(option_Filter);
    	    options.addOption(option_MaxLength);
    	    options.addOption(option_minPercentIdent);
    	    options.addOption(option_minComplexity);
    	    options.addOption(option_Alignment);
    	    
    	    options.addOption(option_Help);
    	    options.addOption(option_Path);
    	    options.addOption(option_Verbose);
    	    options.addOption(option_Reads);
    	    //options.addOption(option_Crawl);
    	    //options.addOption(option_List);
    	    options.addOption(option_MeganSummaries);
    	    options.addOption(option_DeStackOff);
    	    options.addOption(option_DeDupOff);
    	    options.addOption(option_Downsampling);

    	    try //evaluate commandline
    	    {
    	        commandLine = parser.parse(options, parameters);

    	        if (commandLine.hasOption("input"))//evaluate input directorty
    	        {
    	        	log.log(Level.INFO,"Input Set to: ");
    	            for(String arg :commandLine.getOptionValues("input")){
    	            	try{
    	            		File inFile = null;
    	            		if(new File(arg).getParent()!=null){//get Path of inFile either by getting the canonical Path from the parent
    	            			inFile = new File(new File(arg).getParentFile().getCanonicalPath()+"/"+ new File(arg).getName());
    	            		}else {// or by pathching together the Path
    	            			inFile = new File(System.getProperty("user.dir")+"/"+arg);
    	            		}
    	            		if(inFile.isDirectory()){ // if the file is an directory
    	            			 log.info(arg);
    	            			for(String name : inFile.list())//if file ends with RMA6 or is as a soft link at to files
    	            				if(name.endsWith("rma6")|| Files.isSymbolicLink(new File(inFile.getPath()+"/" + name).toPath()))
    	            				 fileNames.add(inFile.getPath()+"/" + name);
    	            		}else if(inFile.isFile()){// is File
    	            			if(arg.endsWith("rma6")||Files.isSymbolicLink(new File(inFile.getPath()).toPath())){ // test if file can be read as if text file that contains names
    	            				log.info(inFile.getPath());
    	            				fileNames.add(inFile.getPath());
    	            			}else{
    	            			readFileList(inFile);}
    	            		}
    	            	}catch(IOException io){
    	            		warning.log(Level.SEVERE,"Can't open File", io);
    	            	}	
    	            }  
    	        }
    	        
    	        if (commandLine.hasOption("output"))//set output directorty
    	        {
    	        	log.log(Level.INFO,"Output Directory set to: "+commandLine.getOptionValue("output"));
    	        	try{
    	        		String path = commandLine.getOptionValue("output");
    	        		if(path.contains("\n")||path.contains("$")||path.contains("\'")||path.contains("=")|| path.contains("\"")){
    	        			warning.log(Level.SEVERE,"Illegal character used in path for out dir");
    	        			System.exit(1);
    	        		}
    	        		File f = new File(path);
    	        		f.isDirectory();
    	        		 outDir = f.getCanonicalPath()+"/";
    	        		}catch(IOException io){
    	        			warning.log(Level.SEVERE,"no valid outdir", io);
    	        		}
    	        }
    	        
    	        if (commandLine.hasOption("taxa"))//evaluate node list
    	        {	 taxas = Taxas.USER;
    	            for(String tax : commandLine.getOptionValues("taxa")){
    	            	log.info("Taxon File: ");
    	            	log.info(tax);
    	            	try{
    	            		File f = new File(tax);
    	     				if(f.getCanonicalFile().exists()){
    	     					readTaxList(f);
    	     				}else{
							   log.info("Added Taxon: ");
							   log.info(tax +" to analysis");
							   taxNames.add(tax.replace('_', ' ')); 
    	     				}
    	            	} catch (IOException e) {
    	            	System.exit(1);
						warning.log(Level.WARNING,"IOException taxa list cannot be resolved",e);
					}
    	           }	   
    			}
    					        
    	        
    	        if (commandLine.hasOption("threads"))//set set number of threads
    	        {	
    	        	log.info("Trying to use "+commandLine.getOptionValue("threads")+" threads");
    	             numThreads=Integer.parseInt(commandLine.getOptionValue("threads"));
    	            if( numThreads > Runtime.getRuntime().availableProcessors()){
    	    			 numThreads = Runtime.getRuntime().availableProcessors(); // enforce that not more processors than available to jvm should be used 
    	    			log.warning("The Number of Threads higher than than Number available to System");
    	            }
    	    		log.log(Level.INFO,"Using " + numThreads +" threads");
    	        }

    	        if (commandLine.hasOption("filter"))// set filters to trigger specific behaviour 
    	        {
    	            if(Pattern.compile(Pattern.quote("default"), Pattern.CASE_INSENSITIVE).matcher(commandLine.getOptionValue("filter")).find()){
	    	            	 behave = Filter.NON;
	    	            	log.log(Level.INFO,"Custom Behaviour set to: "+commandLine.getOptionValue("filter"));
    	            }
    	            else if(Pattern.compile(Pattern.quote("ancient"), Pattern.CASE_INSENSITIVE).matcher(commandLine.getOptionValue("filter")).find()){
    	            		 behave = Filter.ANCIENT;
    	            		log.log(Level.INFO,"Custom Behaviour set to: "+commandLine.getOptionValue("filter"));
    	            	}else if(Pattern.compile(Pattern.quote("scan"), Pattern.CASE_INSENSITIVE).matcher(commandLine.getOptionValue("filter")).find()){
	    	            	 behave = Filter.SCAN;
	    	            	log.log(Level.INFO,"Custom Behaviour set to: "+commandLine.getOptionValue("filter"));
    	            }else if(Pattern.compile(Pattern.quote("crawl"), Pattern.CASE_INSENSITIVE).matcher(commandLine.getOptionValue("filter")).find()){
	        	        	 behave = Filter.CRAWL;
	        	        	log.log(Level.INFO,"Custom Behaviour set to: "+commandLine.getOptionValue("filter"));
        	        }else if(Pattern.compile(Pattern.quote("def_anc"), Pattern.CASE_INSENSITIVE).matcher(commandLine.getOptionValue("filter")).find()){
	    	            	 behave = Filter.NON_ANCIENT;
	    	            	log.log(Level.INFO,"Custom Behaviour set to: "+commandLine.getOptionValue("filter"));
    	            }else if(Pattern.compile(Pattern.quote("full"), Pattern.CASE_INSENSITIVE).matcher(commandLine.getOptionValue("filter")).find()){//TODO remove in future
	    	            	 behave = Filter.NON_ANCIENT;
	    	            	log.log(Level.INFO,"Custom Behaviour set to: "+commandLine.getOptionValue("filter"));
    	            }else if(Pattern.compile(Pattern.quote("srna"), Pattern.CASE_INSENSITIVE).matcher(commandLine.getOptionValue("filter")).find()){//TODO remove in future
	    	            	 behave = Filter.SRNA;
	    	            	log.log(Level.INFO,"Custom Behaviour set to: "+commandLine.getOptionValue("filter"));
    	            }
    	        }

    	        if (commandLine.hasOption("top"))
    	        {
    	        	log.log(Level.INFO,"Top Percent value set to: "+commandLine.getOptionValue("top"));
    	             topPercent = Double.parseDouble(commandLine.getOptionValue("top"));
    	            if(topPercent>1){
    	            	warning.log(Level.SEVERE, "Top percent values higher >1 not accepted please use 0.01 to use the highest scoring percent of alingments");
    	            	System.exit(1);
    	            }
    	            if(topPercent==1){
    	            	log.log(Level.WARNING, "Using top percent of 1 will use all alignmets per read please consider se 0.01 to use the highest scoring percent of alingments");
    	            }
    	        }
    	        
    	        if (commandLine.hasOption("maxReadLength"))
    	        {
    	        	log.log(Level.INFO,"maximum Read Length set to: " + commandLine.getOptionValue("maxReadLength"));
    	             maxLength = Integer.parseInt(commandLine.getOptionValue("maxReadLength"));
    	        }
    	        
    	        if(commandLine.hasOption("minPI")){
    	        	double minP = Double.parseDouble(commandLine.getOptionValue("minPI"));
    	        	if(minP < 1.0){
    	        		log.log(Level.WARNING, "Minimum pIdent "+ minP+" smaller than 1.0");
    	        		minP *= 100;
    	        		log.log(Level.WARNING, "Set Minimum pIdent "+ minP);
    	        		 minPIdent = minP;
    	        	}else{
    	        		 minPIdent = minP;
    	        		log.log(Level.INFO, "Minimum pIdent set to "+ minPIdent);
    	        	}
    	        }
    	        
    	        if(commandLine.hasOption("resources")){
    	        	File f = new File(commandLine.getOptionValue("resources"));
    	        	if(f.isDirectory()){
    	        		 tree_Path = commandLine.getOptionValue("resources");
    	        	}
    	        }
    	        if(commandLine.hasOption("minComp")){
    	        	 minComplexity = Double.parseDouble(commandLine.getOptionValue("minComp"));
    	        	log.log(Level.INFO, "Minimum Complexity set to " + minComplexity);
    	        }
    	        if(commandLine.hasOption('v')){
    	        	log.log(Level.INFO, "Verbose");
    	        	 verbose = true;
    	        }
    	        
    	        if((commandLine.hasOption("reads") && behave != Filter.SCAN)){
    	        	log.log(Level.INFO, "Retrieve filtered Reads");
    	        	 reads = true;
    	        }
    	        if(commandLine.hasOption("list")){
    	        	String list = commandLine.getOptionValue("list");
    	        	String line = "";
    	        	if(Pattern.compile(Pattern.quote("default"), Pattern.CASE_INSENSITIVE).matcher(list).find()){
    	        		log.log(Level.INFO, "use default list" + minComplexity);
    	        		line = "myPointer";
    	        	}// use other lists here
    	        	try{
	            		File f = new File(line);
	     				if(f.getCanonicalFile().exists()){
	     					readTaxList(f);
	     				}
    	        	} catch(IOException io)	{
    	        		warning.log(Level.WARNING, line+" not exist", io);
    	        	}
    	        }
    	        if(commandLine.hasOption("matches")&& behave != Filter.SCAN){
    	        	log.log(Level.INFO, "retrieve Alignments");
    	        	alignment = true;
    	        }
    	        if(commandLine.hasOption("meganSummary")){//request megan summary files
    	        	wantMeganSummaries = true;
    	        }
    	        if(commandLine.hasOption("destackingOff")){
    	        	 destackOff = true;
    	        }
    	        if(commandLine.hasOption("dupRemOff")){
    	        	 deDupOff = true;
    	        }
    	        if(commandLine.hasOption("downSampOff")){
    	        	 downsampling = false;
    	        }
    	        if(commandLine.hasOption("h")){////help
    	        	String header = "MaltExtract beta version 1.3";
    	    	    String footer = "In case you encounter an error drop an email with an useful description to huebler@shh.mpg.de ";
    	    	    HelpFormatter formatter = new HelpFormatter();
    	    	    formatter.setWidth(500);
    	    	    formatter.printHelp("RMAExtractor", header, options, footer, true);   
    	    	    System.exit(0);
    	        }
    	        if(!commandLine.hasOption("h")  && (!commandLine.hasOption("output") || !commandLine.hasOption("input"))){// use if output is missing
    	        	warning.log(Level.SEVERE,"Please, specifiy input files or input directories and output directory");
    	        	System.exit(1);
    	        }
    	        String[] remainder = commandLine.getArgs();// check if unrecognized parameters are missing 
    	        if(remainder.length != 0){
    	           warning.log(Level.SEVERE,"Remaining arguments: ");
    	            for (String argument : remainder){
    	              warning.log(Level.WARNING,argument);
    	            }
    	          }
    	        if(!commandLine.hasOption("taxa") && behave != Filter.SCAN){//stop if taxa list is missing
    	        	warning.log(Level.SEVERE,"No target species provided for filtering");
    	        	System.exit(1);
    	        }//stop if parameters are missing
    	        if(commandLine.hasOption("taxa") && commandLine.hasOption("list")){
    	        	warning.log(Level.SEVERE,"Use either list for build in lists or use taxa not both");
    	        	System.exit(1);
    	        }
    	    }
    	    catch (ParseException exception){
    	    	warning.log(Level.SEVERE, "Parse error: " + exception);
    	    }
    	    
    }
  
}