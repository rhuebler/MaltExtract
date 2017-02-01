package RMA6OutputProcessor;

import java.io.File;
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
import RMA6TaxonProcessor.NodeProcessor;
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
	HashMap<Integer,Integer> defaultSum;
	HashMap<Integer,Integer> ancientSum;
	HashMap<Integer,Integer> nonDuplicateSum;
	HashMap<Integer,Integer> ancientNonDuplicateSum;
	private boolean alignment = false;
	private boolean reads = false;
	// constructor
	public RMA6OutputProcessor(String fileName, String outDir, NCBI_MapReader mapReader,
			Logger warning, Filter behave, boolean alignment, boolean reads) {
		this.fileName = fileName;
		this.outDir = outDir;
		this.mapReader = mapReader;
		this.warning = warning;
		this.behave = behave;
		this.alignment = alignment;
		this.reads =  reads;
	}
	private void setSumLine(HashMap<Integer,Integer> list)
	{	this.defaultSum = list;
	}
	public HashMap<Integer,Integer> getSumLine(){
		return this.defaultSum;
	}
	private void setAncientLine(HashMap<Integer,Integer> list){
		this.ancientSum = list;
	}
	public HashMap<Integer,Integer> getAncientLine(){
		return this.ancientSum;
	}
	private void setNonDuplicateLine(HashMap<Integer,Integer> list){
		this.nonDuplicateSum = list;
	}
	public HashMap<Integer,Integer> getNonDuplicateLine(){
		return this.nonDuplicateSum;
	}
	private void setAncientNonDuplicateLine(HashMap<Integer,Integer> list){
		this.ancientNonDuplicateSum = list;
	}
	public HashMap<Integer,Integer> getAncientNonDuplicateLine(){
		return this.ancientNonDuplicateSum;
	}
	private void writeFilter(ArrayList<String> summary, String outDir){
		String header = "Node\tNumberOfUnfilteredReads\tNumberOfFiltedReads\tNumberOfUnfiltedMatches\tnumberOfMatches";
		summary.sort(null);
		summary.add(0,header);
	try{
		Path file = Paths.get(outDir);
		Files.write(file, summary, Charset.forName("UTF-8"));
	}catch(IOException io){
		warning.log(Level.SEVERE, "File writing exception", io);
	}catch(NullPointerException np){
		warning.log(Level.SEVERE,"Cannot write file"+fileName, np);
	}
	}
	private void writeMisMap(ArrayList<String> summary, String outDir){
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
			Path file = Paths.get(outDir);
			Files.write(file, summary, Charset.forName("UTF-8"));
		}catch(IOException io){
			warning.log(Level.SEVERE, "File writing exception", io);
		}catch(NullPointerException np){
			warning.log(Level.SEVERE,"Cannot write file"+fileName, np);
		}
	 
	}
	private void writeBlastHits(int taxID, ArrayList<String> summary, String outDir){
		try{
			String name;
			if(mapReader.getNcbiIdToNameMap().get(taxID) != null){
				name = mapReader.getNcbiIdToNameMap().get(taxID).replace(' ', '_');
				name = name.replace("/", "_");
				name = name.replace("\\", "_");
				name = name.replace("#", "_");
			}else{
				name = "unassingned_name";
			}	
		Path file = Paths.get(outDir+name+".txt");
		
		Files.write(file, summary, Charset.forName("UTF-8"));
		}catch(IOException io){
			warning.log(Level.SEVERE,"Cannot write file"+ fileName, io);
		}catch(NullPointerException np){
			warning.log(Level.SEVERE,"Cannot write file"+fileName, np);
		}
	}
	private void writeReads(int taxID, ArrayList<String> summary, String outDir){
		try{//TODO look more into this here 
			String name;
			if(mapReader.getNcbiIdToNameMap().get(taxID) != null){
				name = mapReader.getNcbiIdToNameMap().get(taxID).replace(' ', '_');
				name = name.replace("/", "_");
				name = name.replace("\\", "_");
				name = name.replace("#", "_");
			}else{
				name = "unassingned_name";
			}	
		Path file = Paths.get(outDir+name+"done"+".fasta");
		Files.write(file, summary, Charset.forName("UTF-8"));
		}catch(IOException io){
			warning.log(Level.SEVERE,"Cannot write file"+ fileName, io);
		}catch(NullPointerException np){
			warning.log(Level.SEVERE,"Cannot write file"+fileName, np);
		}
	}
	private void writeCoverageHistogram(List<String> summary, String outDir){
		try{
			summary.sort(null);
			String header = "Taxon\tReference\t0\t1\t2\t3\t4\t5\t6\t7\t8\t9\t10\thigher";
			summary.add(0, header);
			Path file = Paths.get(outDir);
			Files.write(file, summary, Charset.forName("UTF-8"));
		}catch(IOException io){
			warning.log(Level.SEVERE,"Cannot write file "+ fileName, io);
		}catch(NullPointerException np){
			warning.log(Level.SEVERE,"Cannot write file"+fileName, np);
		}
	}
	private void writeReadDist(List<String> summary, String outDir){
			//problem is null lines in file
			try{
				summary.sort(null);
				String header = "Taxon\tReference\tuniquePerReference\tnonStacked\tnonDuplicatesonReference\tTotalReadsOnReference\tReferenceLength";
				summary.add(0, header);
				Path file = Paths.get(outDir);
				Files.write(file, summary, Charset.forName("UTF-8"));
			}catch(IOException io){
				warning.log(Level.SEVERE,"Cannot write file " + outDir, io);
			}catch(NullPointerException np){
				warning.log(Level.SEVERE,"Cannot write file "+ outDir, np);
			}
	}
	private void writeEditDistance(List<String> histo, String outDir){
		try{
			String header = "Node\t0\t1\t2\t3\t4\t5\thigher";
			histo.sort(null);
			histo.add(0,header);
			Path file = Paths.get(outDir);
			Files.write(file, histo, Charset.forName("UTF-8"));
		}catch(IOException io){
			warning.log(Level.SEVERE,"Cannot write file"+fileName, io);
		}catch(NullPointerException np){
			warning.log(Level.SEVERE,"Cannot write file"+fileName, np);
		}
	}
	private void writePercentIdentity(List<String> histo, String outDir){
		try{
			String header = "Node\t80\t85\t90\t95\t100";
			histo.sort(null);
			histo.add(0,header);
			Path file = Paths.get(outDir);
			Files.write(file, histo, Charset.forName("UTF-8"));
		}catch(IOException io ){
			warning.log(Level.SEVERE,"Cannot write file", io);
		}catch(NullPointerException np){
			warning.log(Level.SEVERE,"Cannot write file"+fileName, np);
		}
	} 
	private void writeReadLengthDistribution(List<String> histo, String outDir){
		try{
			String header = "Node\tMean\tGeometricMean\tMedian\tStandardDev";
			histo.sort(null);
			histo.add(0,header);
			Path file = Paths.get(outDir);
			Files.write(file, histo, Charset.forName("UTF-8"));
		}catch(IOException io){
			warning.log(Level.SEVERE,"Cannot write file", io);
		}
	}
	private void prepareOutput(HashMap<Integer,Future<NodeProcessor>> results, Filter switcher){
		
		HashMap<Integer,Integer> overallSum = new HashMap<Integer,Integer>();
		ArrayList<String> editDistance = new ArrayList<String>();
		ArrayList<String> percentIdentity = new ArrayList<String>();
		ArrayList<String> readDistribution = new ArrayList<String>();
		ArrayList<String> misMatches = new ArrayList<String>();
		ArrayList<String> coverageHistogram = new ArrayList<String>();
		ArrayList<String> readLengthHistogram = new ArrayList<String>();
		ArrayList<String> filterInformation = new ArrayList<String>();
		for(int id : results.keySet()){
			RMA6TaxonProcessor taxProcessor;
			try {
				if(switcher == Filter.NON){
					taxProcessor = results.get(id).get().getDefault();
				}else if(switcher == Filter.ANCIENT){
					taxProcessor = results.get(id).get().getAncient();
				}else if(switcher == Filter.NONDUPLICATES){
					taxProcessor = results.get(id).get().getNonDuplicate();
				}else {
					taxProcessor = results.get(id).get().getAncientNonDuplicate();
				}
				overallSum.put(id,taxProcessor.getNumberOfReads());
				readDistribution.add(taxProcessor.getReadDistribution());	
				editDistance.add(taxProcessor.getEditDistanceHistogram());
				percentIdentity.add(taxProcessor.getPercentIdentityHistogram());
				coverageHistogram.add(taxProcessor.getCoverageLine());
				misMatches.add(taxProcessor.getDamageLine());
				readLengthHistogram.add(taxProcessor.getReadLengthStatistics());
				filterInformation.add(taxProcessor.getFilterLine());
				if(switcher == Filter.ALL && alignment )
					writeBlastHits(id,taxProcessor.getAlignments(),outDir+"/ancientNonDuplicates/alignments/"+fileName+"/");
				if(switcher == Filter.ANCIENT && alignment)
					writeBlastHits(id,taxProcessor.getAlignments(),outDir+"/ancient/alignments/"+fileName+"/");
				if(switcher == Filter.NON && alignment)
					writeBlastHits(id,taxProcessor.getAlignments(),outDir+"/default/alignments/"+fileName+"/");
				if(switcher == Filter.NONDUPLICATES && alignment)
					writeBlastHits(id,taxProcessor.getAlignments(),outDir+"/nonDuplicates/alignments/"+fileName+"/");
				
				if(switcher == Filter.ALL && reads)
					writeReads(id,taxProcessor.getReads(),outDir+"/ancientNonDuplicates/reads/"+fileName+"/");
				if(switcher == Filter.ANCIENT && reads)
					writeReads(id,taxProcessor.getReads(),outDir+"/ancient/reads/"+fileName+"/");
				if(switcher == Filter.NON && reads)
					writeReads(id,taxProcessor.getReads(),outDir+"/default/reads/"+fileName+"/");
				if(switcher == Filter.NONDUPLICATES && reads)
					writeReads(id,taxProcessor.getReads(),outDir+"/nonDuplicates/reads/"+fileName+"/");
				
			} catch (InterruptedException | ExecutionException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}	
			String dir = "";
			if(switcher == Filter.NON){
				dir = outDir+"/default/";
				setSumLine(overallSum);
			}else if(switcher == Filter.ANCIENT){
				setAncientLine(overallSum);
				dir = outDir+"/ancient/";
			}else if(switcher == Filter.NONDUPLICATES){
				dir = outDir+"/nonDuplicates/";
				setNonDuplicateLine(overallSum);
			}else if(switcher == Filter.ALL){
				dir = outDir+"/ancientNonDuplicates/";
				setAncientNonDuplicateLine(overallSum);
			}
			//write output
			writeMisMap(misMatches, dir+"damageMismatch/"+fileName+"_damageMismatch"+".txt");
			writeReadDist(readDistribution, dir+"readDist/"+fileName+"_readDist"+".txt");
			writeEditDistance(editDistance, dir+"editDistance/"+fileName+"_editDistance"+".txt");
			writePercentIdentity(percentIdentity, dir+"percentIdentity/"+fileName+"_percentIdentity"+".txt");
			writeCoverageHistogram(coverageHistogram, dir+"readDist/"+fileName+"_coverageHist"+".txt");	
			writeReadLengthDistribution(readLengthHistogram, dir+"readDist/"+fileName+"_readLengthDist"+".txt");
			writeFilter(filterInformation,dir+"FilterInformation/"+fileName+"_filterTable"+".txt");
	}
	public void process(HashMap<Integer,Future<NodeProcessor>> results){
		if(behave == Filter.NON){
			if(alignment){
				new File(outDir+"/default/"+"/alignments/"+fileName).mkdirs();
			}
			if(reads){
				new File(outDir+"/default/"+"/reads/"+fileName).mkdirs();
			}
			prepareOutput(results,behave);
		}else if(behave == Filter.ANCIENT){
			if(alignment){
				new File(outDir+"/ancient/"+"/alignments/"+fileName).mkdirs();
			}
			if(reads){
				new File(outDir+"/ancient/"+"/reads/"+fileName).mkdirs();
			}
			prepareOutput(results,behave);
		}else if(behave == Filter.NONDUPLICATES){
			if(alignment){
				new File(outDir+"/nonDuplicates/"+"/alignments/"+fileName).mkdirs();
			}
			if(reads){
				new File(outDir+"/nonDuplicates/"+"/reads/"+fileName).mkdirs();
			}
			prepareOutput(results,behave);
		}else if(behave == Filter.ALL){
			if(alignment){
				new File(outDir+"/ancientNonDuplicates/"+"/alignments/"+fileName).mkdirs();
			}
			if(reads){
				new File(outDir+"/ancientNonDuplicates/"+"/reads/"+fileName).mkdirs();	
			}
			prepareOutput(results,behave);
		}else if(behave == Filter.NON_ANCIENT){
			
			if(alignment){
				new File(outDir+"/ancient/"+"/alignments/"+fileName).mkdirs();
				new File(outDir+"/default/"+"/alignments/"+fileName).mkdirs();
			}	
			if(reads){
					new File(outDir+"/default/"+"/reads/"+fileName).mkdirs();
					new File(outDir+"/ancient/"+"/reads/"+fileName).mkdirs();
			}
			prepareOutput(results,Filter.NON);
			prepareOutput(results,Filter.ANCIENT);
		}
	}
}
