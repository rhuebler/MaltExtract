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
import behaviour.OutputType;
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
	private void setSumLine(HashMap<Integer,Integer> list)// set Sumline
	{	this.defaultSum = list;
	}
	public HashMap<Integer,Integer> getSumLine(){//get SumLine
		return this.defaultSum;
	}
	private void setAncientLine(HashMap<Integer,Integer> list){// get ancient sum Sumline
		this.ancientSum = list;
	}
	public HashMap<Integer,Integer> getAncientLine(){ //get Ancient Line
		return this.ancientSum;
	}
	private void setNonDuplicateLine(HashMap<Integer,Integer> list){// non duplicate line
		this.nonDuplicateSum = list;
	}
	public HashMap<Integer,Integer> getNonDuplicateLine(){// get nonduplicate line
		return this.nonDuplicateSum;
	}
	private void setAncientNonDuplicateLine(HashMap<Integer,Integer> list){// set ancient non duplicate line
		this.ancientNonDuplicateSum = list;
	}
	public HashMap<Integer,Integer> getAncientNonDuplicateLine(){// get ancient non duplicate line
		return this.ancientNonDuplicateSum;
	}
	String getName(int taxID){
		String name;
		if(mapReader.getNcbiIdToNameMap().get(taxID) != null){
			name = mapReader.getNcbiIdToNameMap().get(taxID).replace(' ', '_');
			name = name.replace("/", "_");
			name = name.replace("\\", "_");
			name = name.replace("#", "_");
		}else{
			name = "unassingned_name";
		}	
		return name;
	}
	private void writeOutput(List<String> histo, String dir,OutputType type,int taxID){
		try{
			String  header ="Taxon";
			String outDir = dir;
			switch (type) {
			case ADDIDTIONALENTRIES: header = "TargetNode\t1.0\t0.9\t0.8\t0.7\t0.6\t0.5\t0.4\t0.3\t0.2\t0.1";
									outDir += "readDist/"+fileName+"_additionalNodeEntries"+".txt";
									break;
			case ALIGNMENTS:		outDir += getName(taxID)+".fasta";
									break;
			case ALIGNMENTDISTRIBUTION: header = "Taxon\tReference\tuniquePerReference\tnonStacked\tnonDuplicatesonReference\tTotalAlignmentsOnReference\tReferenceLength\tAverageCoverage";
										outDir+= "readDist/"+fileName+"_alignmentDist"+".txt";
									break;
			case COVERAGEHISTOGRAM: header = "Taxon\tReference\t0\t1\t2\t3\t4\t5\t6\t7\t8\t9\t10\thigher";
									outDir += "coverage/"+fileName+"_coverageHist"+".txt";
									break;
			case DAMAGE: header = "Node";
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
						 	outDir += "damageMismatch/"+fileName+"_damageMismatch"+".txt";
						 			break;
			case EDITDISTANCE: header="Taxon\tReference\t0\t1\t2\t3\t4\t5\t6\t7\t8\t9\t10\thigher";
								outDir += "editDistance/"+fileName+"_editDistance"+".txt";
									break;
			case FILTER: header= "Node\tNumberOfUnfilteredReads\tNumberOfFilteredReads\tNumberOfUnfilteredAlignments\tnumberOfAlignments\tturnedOn?";
						 outDir += "FilterInformation/"+fileName+"_filterTable"+".txt";			
									break;
			case READS:	outDir += getName(taxID)+".fasta";
						break;
			case READLENGTH_DIST: for(int i = 25;i<=200;i+=5){
									header+="\t"+i+"bp";
									}
									outDir+="readDist/"+fileName+"_readLengthDist"+".txt";
									break;
			case READLENGTH_STATISTICS:	header = "Taxon\tMean\tGeometricMean\tMedian\tStandardDev";
									outDir+="readDist/"+fileName+"_readLengthStat"+".txt";
									break;
			
			case PERCENTIDENTITY:	header = "Taxon\t80\t85\t90\t95\t100";
									outDir += "percentIdentity/"+fileName+"_percentIdentity"+".txt";
									break;
			case POS_COVERED:	header = "Taxon\tReference\tAverageCoverage\tCoverge_StandardDeviation\tpercCoveredHigher1\tpercCoveredHigher2\tpercCoveredHigher3\tpercCoveredHigher4\tpercCoveredHigher5";
								outDir += "coverage/"+fileName+"_postionsCovered"+".txt";
									break;
			default:
				warning.log(Level.SEVERE,"No output specified");
				break;
			};
			histo.sort(null);
			histo.add(0,header);
			Path file = Paths.get(outDir);
			Files.write(file, histo, Charset.forName("UTF-8"));
		}catch(IOException io){
			warning.log(Level.SEVERE,"Cannot write file", io);
		}
	}
	private void prepareOutput(HashMap<Integer, Future<NodeProcessor>> hashMap, Filter switcher){// Gather information from node list and write ouput by preparing the information
		HashMap<Integer,Integer> overallSum = new HashMap<Integer,Integer>();
		
		ArrayList<String> additionalEntries = new ArrayList<String>();
		ArrayList<String> coverageHistogram = new ArrayList<String>();
		ArrayList<String> coveragePositions = new ArrayList<String>();
		ArrayList<String> editDistance = new ArrayList<String>();
		ArrayList<String> misMatches = new ArrayList<String>();
		ArrayList<String> percentIdentity = new ArrayList<String>();
		ArrayList<String> readDistribution = new ArrayList<String>();
		ArrayList<String> readLengthHistogram = new ArrayList<String>();
		ArrayList<String> readLengthStatistics = new ArrayList<String>();
		ArrayList<String> filterInformation = new ArrayList<String>();
				
		String dir = "";
		if(switcher == Filter.NON){
			dir = outDir+"default/";
			setSumLine(overallSum);
		}else if(switcher == Filter.ANCIENT){
			setAncientLine(overallSum);
			dir = outDir+"ancient/";
		}else if(switcher == Filter.NONDUPLICATES){
			dir = outDir+"nonDuplicates/";
			setNonDuplicateLine(overallSum);
		}else if(switcher == Filter.ALL){
			dir = outDir+"ancientNonDuplicates/";
			setAncientNonDuplicateLine(overallSum);
		}
		
		for(int id : hashMap.keySet()){
			try{
				NodeProcessor nodeProcessor = hashMap.get(id).get();
				RMA6TaxonProcessor taxProcessor;
					if(switcher == Filter.NON){
						taxProcessor = nodeProcessor.getDefault();
					}else if(switcher == Filter.ANCIENT){
						taxProcessor = nodeProcessor.getAncient();
					}else{
						taxProcessor=null;
						warning.log(Level.SEVERE, "Filter no longer supported");
						System.exit(1);
					}
					overallSum.put(id,taxProcessor.getNumberOfReads());
					additionalEntries.add(taxProcessor.getAdditionalEntries());
					coverageHistogram.add(taxProcessor.getCoverageLine());
					coveragePositions.add(taxProcessor.getCoveragePositions());
					editDistance.add(taxProcessor.getEditDistanceHistogram());
					filterInformation.add(taxProcessor.getFilterLine());
					misMatches.add(taxProcessor.getDamageLine());
					percentIdentity.add(taxProcessor.getPercentIdentityHistogram());
					readDistribution.add(taxProcessor.getReadDistribution());	
					readLengthHistogram.add(taxProcessor.getReadLengthDistribution());
					readLengthStatistics.add(taxProcessor.getReadLengthStatistics());	
					if(alignment )
						writeOutput(taxProcessor.getAlignments(),dir+"/alignments/"+fileName+"/",OutputType.ALIGNMENTS, id);
					if(reads)
						writeOutput(taxProcessor.getReads(),dir+"/reads/"+fileName+"/",OutputType.READS, id);
			}catch(InterruptedException | ExecutionException  e){
				warning.log(Level.SEVERE,"Error", e);
			}
		}
			
			//write output
			writeOutput(additionalEntries, dir,OutputType.ADDIDTIONALENTRIES, 0);
			writeOutput(readDistribution, dir,OutputType.ALIGNMENTDISTRIBUTION, 0);
			writeOutput(coveragePositions, dir,OutputType.POS_COVERED, 0);
			writeOutput(coverageHistogram, dir,OutputType.COVERAGEHISTOGRAM, 0);
			writeOutput(misMatches, dir,OutputType.DAMAGE, 0);
			writeOutput(editDistance, dir,OutputType.EDITDISTANCE, 0);
			writeOutput(filterInformation,dir,OutputType.FILTER, 0);
			writeOutput(percentIdentity, dir,OutputType.PERCENTIDENTITY, 0);
			writeOutput(readLengthHistogram, dir,OutputType.READLENGTH_DIST, 0);
			writeOutput(readLengthStatistics, dir, OutputType.READLENGTH_STATISTICS, 0);
			hashMap.clear();
	}
	public void process(HashMap<Integer, Future<NodeProcessor>> hashMap){ // process NodeProcessorts depending on files 
		if(behave == Filter.NON){// retrieve the filter that was used and process accordingly
			// do not use ancient filter
			if(alignment){
				new File(outDir+"/default/"+"/alignments/"+fileName).mkdirs();
			}
			if(reads){
				new File(outDir+"/default/"+"/reads/"+fileName).mkdirs();
			}
			prepareOutput(hashMap,behave);
		}else if(behave == Filter.ANCIENT){// just use ancient filter
			if(alignment){
				new File(outDir+"/ancient/"+"/alignments/"+fileName).mkdirs();
			}
			if(reads){
				new File(outDir+"/ancient/"+"/reads/"+fileName).mkdirs();
			}
			prepareOutput(hashMap,behave);
		}else if(behave == Filter.NON_ANCIENT){// more or less the default case
			
			if(alignment){
				new File(outDir+"/ancient/"+"/alignments/"+fileName).mkdirs();
				new File(outDir+"/default/"+"/alignments/"+fileName).mkdirs();
			}	
			if(reads){
					new File(outDir+"/default/"+"/reads/"+fileName).mkdirs();
					new File(outDir+"/ancient/"+"/reads/"+fileName).mkdirs();
			}
			prepareOutput(hashMap,Filter.NON);
			prepareOutput(hashMap,Filter.ANCIENT);
		}
	}
}
