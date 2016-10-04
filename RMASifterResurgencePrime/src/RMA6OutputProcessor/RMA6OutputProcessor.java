package RMA6OutputProcessor;

import java.io.IOException;
import java.nio.charset.Charset;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.Future;
import java.util.logging.Level;
import java.util.logging.Logger;

import NCBI_MapReader.NCBI_MapReader;
import RMA6TaxonProcessor.RMA6TaxonProcessor;
import behaviour.Filter;

public class RMA6OutputProcessor {
	private String fileName;
	private String outDir;
	private NCBI_MapReader mapReader;
	private Logger warning;
	private Filter behave;
	HashMap<Integer,Integer> overallSum;
	private boolean reads = false;
	// constructor
	public RMA6OutputProcessor(String fileName, String outDir, NCBI_MapReader mapReader,
			Logger warning, Filter behave, boolean reads) {
		this.fileName = fileName;
		this.outDir = outDir;
		this.mapReader = mapReader;
		this.warning = warning;
		this.behave = behave;
		this.reads = reads;
	}
	private void setSumLine(HashMap<Integer,Integer> list)
	{	this.overallSum = list;
	}
	public HashMap<Integer,Integer> getSumLine(){
		return this.overallSum;
	}
	private void writeMisMap(ArrayList<String> summary){
			String header = "Node";
			for(int i = 0; i < 20; i++){
				if(i<10)
					header+="\t"+"C>T_"+(i+1);
				else
					header+="\t"+"G>A_"+(i+1);
			}
			header+="\tconsidered_Matches";
			summary.sort(null);
			summary.add(0,header);
		try{
			Path file = Paths.get(outDir+"damageMismatch/"+fileName+".txt");
			Files.write(file, summary, Charset.forName("UTF-8"));
		}catch(IOException io){
			warning.log(Level.SEVERE, "File writing exception", io);
		}
	 
	}
	private void writeBlastHits(int taxID, ArrayList<String> summary){
		try{
			String name;
			if(mapReader.getNcbiIdToNameMap().get(taxID) != null)
				name = mapReader.getNcbiIdToNameMap().get(taxID).replace(' ', '_');
			else
				name = "unassingned_name";
		Path file = Paths.get(outDir+"reads/"+fileName+"/"+name+".txt");
		Files.write(file, summary, Charset.forName("UTF-8"));
		}catch(IOException io){
			warning.log(Level.SEVERE,"Cannot write file", io);
		}
	}
	private void writeReadDist(List<String> summary){
		try{
			summary.sort(null);
			String header = "Taxon\tReference\tMeanReadDistance\tMedianReadDistance\tVarianceReadDistance\tStandardDeviationReadDistance\tuniquePerReference\tnonDuplicatesonReference\tTotalReadsOnReference\tReferenceLength";
			summary.add(0, header);
			Path file = Paths.get(outDir+"/readDist/"+fileName+"_readDist"+".txt");
			Files.write(file, summary, Charset.forName("UTF-8"));
		}catch(IOException io){
			warning.log(Level.SEVERE,"Cannot write file", io);
		}
		//this.readDist = null; //delete data to save space after were done potentially the gc should take of it
	}
	private void writeEditDistance(List<String> histo){
		try{
			String header = "Node\t0\t1\t2\t3\t4\t5\thigher";
			histo.sort(null);
			histo.add(0,header);
			Path file = Paths.get(outDir+"/editDistance/"+fileName+"_editDistance"+".txt");
			Files.write(file, histo, Charset.forName("UTF-8"));
		}catch(IOException io){
			warning.log(Level.SEVERE,"Cannot write file", io);
		}
	}
	private void writePercentIdentity(List<String> histo){
		try{
			String header = "Node\t80\t85\t90\t95\t100";
			histo.sort(null);
			histo.add(0,header);
			Path file = Paths.get(outDir+"/percentIdentity/"+fileName+"_percentIdentity"+".txt");
			Files.write(file, histo, Charset.forName("UTF-8"));
		}catch(IOException io){
			warning.log(Level.SEVERE,"Cannot write file", io);
		}
	} 
	public void prepareOutput(HashMap<Integer,Future<RMA6TaxonProcessor>> results){
		HashMap<Integer,Integer> overallSum = new HashMap<Integer,Integer>();
		ArrayList<String> editDistance = new ArrayList<String>();
		ArrayList<String> percentIdentity = new ArrayList<String>();
		ArrayList<String> readDistribution = new ArrayList<String>();
		ArrayList<String> misMatches = new ArrayList<String>();
		for(int id : results.keySet()){
			RMA6TaxonProcessor taxProcessor;
			try {
				taxProcessor = results.get(id).get();
				overallSum.put(id,taxProcessor.getNumberOfMatches());
				readDistribution.add(taxProcessor.getReadDistribution());	
				editDistance.add(taxProcessor.getEditDistanceHistogram());
				percentIdentity.add(taxProcessor.getPercentIdentityHistogram());
				HashMap<Integer, Integer> map = taxProcessor.getMisMap();
				String line;
				if(mapReader.getNcbiIdToNameMap().get(id) != null)
					line = mapReader.getNcbiIdToNameMap().get(id).replace(' ', '_');
				else
					line = "unassingned_name";
				for(int i = 0;i < 20; i++){
					if(map.containsKey(i))
						line+="\t"+map.get(i);
					else	
						line+="\t"+0;
				}
				line += "\t"+map.get(20);
				misMatches.add(line);
				if((behave == Filter.ALL && reads )|| (behave == Filter.ANCIENT && reads)){
					writeBlastHits(id,taxProcessor.getReads());
				}
			} catch (InterruptedException | ExecutionException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}	
			setSumLine(overallSum); // set number of assigned Reads to overall file summary
			writeMisMap(misMatches);
			writeReadDist(readDistribution); // RMA6Processor now saves its own output 
			writeEditDistance(editDistance);
			writePercentIdentity(percentIdentity);
	}
}
