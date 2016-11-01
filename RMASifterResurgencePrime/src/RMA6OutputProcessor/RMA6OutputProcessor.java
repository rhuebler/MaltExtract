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
/**
 * This class is used to retrieve all outputs from RMA6TaxonProcessors
 * Writes RMA6Distributions. Coverage Histograms, damage mismatch Histograms 
 * Edit Histograms and Percent Identity Histograms and potentially Reads
 * @author huebler
 *
 */
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
			String header_part2 ="";
			for(int i = 0; i < 20; i++){
				if(i<10){
					header+="\t"+"C>T_"+(i+1);
					header_part2+="\t"+"D>V(11Substitution)_"+(i+1);
				}else{
					header+="\t"+"G>A_"+(i+1);
					header_part2+="\t"+"H>B(11Substitution)_"+(i+1);
				}
			}
			header+=header_part2+"\tconsidered_Matches";
			summary.sort(null);
			summary.add(0,header);
		try{
			Path file = Paths.get(outDir+"damageMismatch/"+fileName+"_damageMismatch"+".txt");
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
	private void writeCoverageHistogram(List<String> summary){
		try{
			summary.sort(null);
			String header = "Taxon\tReference\t0\t1\t2\t3\t4\t5\t6\t7\t8\t9\t10\thigher";
			summary.add(0, header);
			Path file = Paths.get(outDir+"/readDist/"+fileName+"_coverageHist"+".txt");
			Files.write(file, summary, Charset.forName("UTF-8"));
		}catch(IOException io){
			warning.log(Level.SEVERE,"Cannot write file", io);
		}
	}
	private void writeReadDist(List<String> summary){
		try{
			summary.sort(null);
			String header = "Taxon\tReference\tuniquePerReference\tnonStacked\tnonDuplicatesonReference\tTotalReadsOnReference\tReferenceLength";
			summary.add(0, header);
			Path file = Paths.get(outDir+"/readDist/"+fileName+"_readDist"+".txt");
			Files.write(file, summary, Charset.forName("UTF-8"));
		}catch(IOException io){
			warning.log(Level.SEVERE,"Cannot write file", io);
		}
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
		ArrayList<String> coverageHistogram = new ArrayList<String>();
		for(int id : results.keySet()){
			RMA6TaxonProcessor taxProcessor;
			try {
				taxProcessor = results.get(id).get();
				overallSum.put(id,taxProcessor.getNumberOfReads());
				readDistribution.add(taxProcessor.getReadDistribution());	
				editDistance.add(taxProcessor.getEditDistanceHistogram());
				percentIdentity.add(taxProcessor.getPercentIdentityHistogram());
				coverageHistogram.add(taxProcessor.getCoverageLine());
		
				misMatches.add(taxProcessor.getDamageLine());
				if((behave == Filter.ALL && reads )|| (behave == Filter.ANCIENT && reads)
						|| (behave == Filter.NON && reads)){
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
			writeCoverageHistogram(coverageHistogram);
	}
}
