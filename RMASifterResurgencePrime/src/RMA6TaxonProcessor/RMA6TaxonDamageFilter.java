package RMA6TaxonProcessor;

import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;

import NCBI_MapReader.NCBI_MapReader;
import RMAAlignment.Alignment;
import RMAAlignment.CompositionMap;
import jloda.util.ListOfLongs;
import megan.data.IMatchBlock;
import megan.data.IReadBlock;
import megan.data.IReadBlockIterator;
import megan.data.ReadBlockIterator;
import megan.rma6.ClassificationBlockRMA6;
import megan.rma6.RMA6File;
import megan.rma6.ReadBlockGetterRMA6;
/**
 * Extract all information from one Taxon and save the information in the specified slots to be retrieved
 * in RMA6Processor this taxon processor only processes reads that have at least one match that hints at c-> end substitutions 
 * it will be used to replace the non filtering taxon processor within RMA6Processor when filtering for damaged Reads 
 * @author huebler
 *
 */
public class RMA6TaxonDamageFilter extends RMA6TaxonProcessor{

public RMA6TaxonDamageFilter(int id, NCBI_MapReader reader) {
		super(id, reader);
		// TODO Auto-generated constructor stub
	}

@Override
public void process(RMA6File rma6File, String fileName, double topPercent, int maxLength){ 
	DecimalFormat df = new DecimalFormat("#.###");
	ArrayList<String> supplemantary = new ArrayList<String>();
	// use ReadsIterator to get all Reads assigned to MegantaxID and print top percent to file;
	try{final ClassificationBlockRMA6 block = new ClassificationBlockRMA6("Taxonomy");
    final long start = rma6File.getFooterSectionRMA6().getStartClassification("Taxonomy");
    block.read(start, rma6File.getReader());
    final ListOfLongs list = new ListOfLongs();
        if (block.getSum(taxID) > 0) {
            block.readLocations(start, rma6File.getReader(), taxID, list);
        }
	IReadBlockIterator classIt  = new ReadBlockIterator(list, new ReadBlockGetterRMA6(rma6File, true, true, (float) 1.0,(float) 100.00,false,true));
	if(!classIt.hasNext()){ // check if reads are assigned to TaxID if not print to console and skip
		System.err.println("TaxID: " + taxID +  " not assigned in File " + fileName+"\n");
		setReadDistribution(mapReader.getNcbiIdToNameMap().get(taxID).replace(' ', '_')+"\tNA\t0\t0\t0\t0\t0\t0\t0\t0");
		setSupplementary(new ArrayList<String>(Arrays.asList("0\t0\t0\t0\t0\t0\t"+mapReader.getNcbiIdToNameMap().get(taxID).replace(' ', '_'))));// in case of taxID not being supported add empty Line
		setNumberOfMatches(0);	
	}else{
	// only consider topscoring alignments

	//System.out.println("Processing Taxon "+mapReader.getNcbiIdToNameMap().get(taxID)+" in File " + fileName); 
	HashMap<Integer, ArrayList<Alignment>> taxonMap = new HashMap<Integer,ArrayList<Alignment>>();
	int numReads = 0;
		while(classIt.hasNext()){
			IReadBlock current = classIt.next();
			if(current.getReadLength()<= maxLength || maxLength == 0){
				IMatchBlock[] blocks=current.getMatchBlocks();
					int k=0;
				float topScore = blocks[0].getBitScore();
				double pIdent = 0;
				double length = 0; 
				int damage=0;
				for(int i = 0; i< blocks.length;i++){
					if(blocks[i].getBitScore()/topScore <= 1 - topPercent){
						break;}
					
					Alignment al = new Alignment();
					al.processText(blocks[i].getText().split("\n"));
					if(al.getFivePrimeDamage()){
						
						length += al.getMlength();
						pIdent += blocks[i].getPercentIdentity();
						if(!taxonMap.containsKey(blocks[i].getTaxonId())){
							ArrayList<Alignment> entry =new ArrayList<Alignment>();
							entry.add(al);
							taxonMap.put(blocks[i].getTaxonId(), entry);
						}else{
							ArrayList<Alignment> entry = taxonMap.get(blocks[i].getTaxonId());
							entry.add(al);
							taxonMap.put(blocks[i].getTaxonId(),entry);
						}
					
						damage++;
					k++;
					}
				}
				if(damage !=0){
					numReads++;
					supplemantary.add(
						current.getReadName()+"\t"
						+ (int)(length/k)+"\t"
						+ df.format(pIdent/(k)) +"\t"
						+ current.getNumberOfMatches()+"\t"
						+ k +"\t"
						+ damage+'\t'
						+ df.format(getGcContent(current.getReadSequence()))+"\t"
						+ mapReader.getNcbiIdToNameMap().get(taxID).replace(' ', '_'));}
			}// if TODO should I add an else here and what to do 
		}// while
			classIt.close();
			String taxName;
			if(mapReader.getNcbiIdToNameMap().get(taxID) != null)
				taxName = mapReader.getNcbiIdToNameMap().get(taxID).replace(' ', '_');
			else
				taxName = "unassingned name";
			CompositionMap map = new CompositionMap(taxonMap);
			map.process();
			String maxReference;
			if(mapReader.getNcbiIdToNameMap().get(map.getMaxID()) != null)
				maxReference =  mapReader.getNcbiIdToNameMap().get(map.getMaxID()).replace(' ', '_');
			else
				maxReference = "unassinged_reference_name";
			String s = taxName + "\t" 
						+ maxReference;
			for(double d : map.getStatistics())
				s += "\t" + df.format(d);
			setReadDistribution(s);
			setNumberOfMatches(numReads);
			setSupplementary(supplemantary);
		 }//else
		}catch(IOException io){
			io.printStackTrace();
		}
	}// void 
}// class 
