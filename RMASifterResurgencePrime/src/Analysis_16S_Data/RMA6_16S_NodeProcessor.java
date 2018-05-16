package Analysis_16S_Data;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.logging.Logger;

import NCBI_MapReader.NCBI_MapReader;
import RMA6TaxonProcessor.RMA6TaxonProcessor;
import RMAAlignment.Alignment;
import RMAAlignment.CompositionMap;
import behaviour.Filter;
import megan.data.IMatchBlock;

public class RMA6_16S_NodeProcessor extends RMA6TaxonProcessor{

	//constructers and set values
	public RMA6_16S_NodeProcessor(int id ,double pID, NCBI_MapReader reader,
			boolean verbose,Logger log, Logger warning, boolean wantReads,double tp,int mL,boolean wantAls,boolean turnOffDestacking,boolean turnOffDeDuping,Filter behave) {
		super(id ,pID, reader, verbose, log, warning, wantReads, tp, mL, wantAls, turnOffDestacking, turnOffDeDuping, behave);
	
	}
	//process each Matchblock
	public void processMatchBlocks(IMatchBlock[] blocks, String name, int length, String sequence){
		originalNumberOfReads++;
		float topScore = blocks[0].getBitScore();
			for(int i = 0; i< blocks.length;i++){
				if(blocks[i].getBitScore()/topScore < 1-topPercent){
					break;}		
				
				Alignment al = new Alignment();
				al.setText(blocks[i].getText());
				al.setTaxID(blocks[i].getTaxonId());
				al.processText();
				al.setPIdent(blocks[i].getPercentIdentity());
				al.setAcessionNumber(blocks[i].getTextFirstWord());
				al.setSequence(sequence);
				al.setReadLength(length);
				al.setReadName(name);
				if(i==0){
					al.setTopAlignment(true);
				}
				originalNumberOfAlignments++;
				if(minPIdent <= al.getPIdent()){ // check for minPercentIdentity
								//get mismatches
					if(!taxonMap.containsKey(al.getTaxID())){
						HashMap<String,ArrayList<Alignment>> list = new HashMap<String,ArrayList<Alignment>>();
						ArrayList<Alignment> entry =new ArrayList<Alignment>();
						entry.add(al);
						list.put(al.getAccessionNumber(), entry);
						taxonMap.put(al.getTaxID(), list);
					}else{
						HashMap<String,ArrayList<Alignment>> list =  taxonMap.get(al.getTaxID());
						if(!list.containsKey(al.getAccessionNumber())){
							ArrayList<Alignment> entry = new ArrayList<Alignment>();
							entry.add(al);
							list.put(al.getAccessionNumber(), entry);
						}else{
							ArrayList<Alignment> entry = list.get(al.getAccessionNumber());
							entry.add(al);
							list.replace(al.getAccessionNumber(),entry);
						}
						
						taxonMap.replace(al.getTaxID(),list);
					}
				}
			}
	}		
	public void process(){ 
		//analyze collected data 
		CompositionMap map = new CompositionMap(taxonMap,turnOffDestacking,turnOffDeDuping, filter);
		map.process();
		map.getNonStacked();
		map.setMapReader(mapReader);
		HashMap<String,Alignment> list = map.getResultsMap();
		for(String key : list.keySet()){
			Alignment al = list.get(key);
			String sequence = "";
			numMatches++;	
			container.processAlignment(al);
			if(wantAlignments){
				alignments.add(al.getReadName());
				alignments.add(al.getText());
			}
	
			lengths.add(al.getReadLength());
			pIdents.add(al.getPIdent());
			distances.add(al.getEditDistance());
			if(wantReads){
				String readName = al.getReadName();
				sequence = al.getSequence();
				if (!readName.startsWith(">"))
					readName = ">"+readName;
				lines.add(readName);
				if (sequence != null) {
					lines.add(sequence);    
				}
			}
		}	
		setOriginalNumberOfAlignments(originalNumberOfAlignments);
		setOriginalNumberOfReads(originalNumberOfReads);
		setDamageLine(container.getDamageLine());
		setNumMatches(numMatches);
		setNumberOfReads(list.keySet().size());
		processCompositionMap(map);
		setEditDistanceHistogram(distances);
		setPercentIdentityHistogram(pIdents);
		setReads(lines);
		setAlignments(alignments);
		calculateReadLengthDistribution();
		map = null;
	}//process
}
