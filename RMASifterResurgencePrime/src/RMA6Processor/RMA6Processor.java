package RMA6Processor;
//import java.io.BufferedOutputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.nio.charset.Charset;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.zip.GZIPOutputStream;

import NCBI_MapReader.NCBI_MapReader;
import NCBI_MapReader.NCBI_TreeReader;
import RMA6TaxonProcessor.RMA6TaxonProcessor;
import megan.rma6.RMA6Connector;
/**
 * Take care of extracting all the information for one RMA6 file and a List of taxons
 * @author huebler
 *
 */
public class RMA6Processor {
	// attributes
	private HashMap<Integer,Integer> overallSum;
	private ArrayList<String> readDist;
	private ArrayList<String> supplement;
	private String outDir;
	private String fileName;
	private String inDir;
	private NCBI_MapReader mapReader;
	private NCBI_TreeReader treeReader;
	private Set<Integer> containedIDs;
	// constructor
	public RMA6Processor(String inDir, String fileName, String outDir, NCBI_MapReader mapReader, NCBI_TreeReader treeReader) {
		this.inDir = inDir;
		this.outDir = outDir;
		this.fileName = fileName;
		this.mapReader = mapReader;
		this.treeReader = treeReader;
	}
	
	//setters
	private void setContainedIDs(Set<Integer> set){
		this.containedIDs = set;
	}
	private void setSumLine(HashMap<Integer,Integer> list)
	{	this.overallSum = list;
	}
	private void setSupplementaryData(ArrayList<String> s){
		String header = "readName\tmatchLength\tPercentIdent\tnumberOfMatches\tconsideredMatches\tdamage\tGC_Content\ttaxName";
		s.add(0, header);
		this.supplement = s ;
	}
	
	private void setReadDistribution(ArrayList<String>r){
		String header = "Taxon\tReference\tMeanReadDistance\tMedianReadDistance\tVarianceReadDistance\tStandardDeviationReadDistance\tuniquePerReference\tnonDuplicatesonReference\tTotalReadsOnReference\tReferenceLength";
		r.add(0,header);
		this.readDist = r;
	}
	// private utility functions
	private void writeReadDist(List<String> summary, String fileName) throws IOException{
		System.out.println("Writing Read Distribution txt File");
		Path file = Paths.get(outDir+fileName+"_readDist"+".txt");
		Files.write(file, summary, Charset.forName("UTF-8"));
		System.out.println("ReadDistribution for " + fileName +" Done!");
	}
	private void writeSupplementary(List<String> supplement, String fileName) throws IOException{ 
		OutputStreamWriter outer = new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(outDir+fileName+"_supplement"+ ".txt.gz")));
		for(String line : supplement)
			outer.write(line);
		outer.close();
		System.out.println("Supplementary for File "+fileName+" done!");
	}
	//getter
	public ArrayList<String> getReadDistribution(){
		return this.readDist;
	}
	public ArrayList<String> getSupplementary(){
		return this.supplement;
	}
	public HashMap<Integer,Integer> getSumLine(){
		return this.overallSum;
	}
	public Set<Integer> getContainedIDs(){
		return this.containedIDs;
	}
	public String getfileName(){
		return this.fileName;
	}
// processing 
public void process(List<Integer>taxIDs, double topPercent) throws IOException{
	RMA6Connector fileCon = new RMA6Connector(inDir+fileName);
	HashMap<Integer,Integer> overallSum = new HashMap<Integer,Integer>();
	ArrayList<String> supplemantary = new ArrayList<String>();
	ArrayList<String> readDistribution = new ArrayList<String>();
	Set<Integer> keys = fileCon.getClassificationBlock("Taxonomy").getKeySet();// get all assigned IDs in a file 
	Set<Integer> idsToProcess = new HashSet<Integer>();
	for(int taxID : taxIDs){
		idsToProcess.add(taxID);
		for(int id :treeReader.getStrains(taxID, keys))
			if(!taxIDs.contains(id))
				idsToProcess.add(id);
	}
	setContainedIDs(idsToProcess);
	for(int id : idsToProcess){
	 RMA6TaxonProcessor taxProcessor = new RMA6TaxonProcessor();// could add new
	 taxProcessor.process(fileCon, id, fileName, mapReader, topPercent);
	 overallSum.put(id,taxProcessor.getNumberOfMatches());
	 readDistribution.add(taxProcessor.getReadDistribution());
	 for(String sup : taxProcessor.getSupplementary())
		 supplemantary.add(sup);
	}//TaxIDs
	setSupplementaryData(supplemantary);	//save supplementary data at read resolution in adequate slot
	setSumLine(overallSum); // set number of assigned Reads to overall file summary
	setReadDistribution(readDistribution);// save ReadDist summary file
	
	writeReadDist(getReadDistribution(),fileName); // RMA6Processor now saves its own output 
	writeSupplementary(getSupplementary(),fileName);
	
 }
}
