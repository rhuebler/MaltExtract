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
	public RMA6_16S_NodeProcessor(int id ,double pID, NCBI_MapReader reader,boolean v,Logger log, Logger warning, boolean reads,double tp,int mL,boolean wantAls, boolean turnOffDestacking,boolean turnOffDeDuping, Filter behave, boolean useAllAlignments, boolean singleStranded) {
		//intialize Nodeanalyzer
		super(id ,pID, reader, v, log, warning, reads, tp, mL, wantAls, turnOffDestacking, turnOffDeDuping, behave, useAllAlignments, singleStranded);
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
				al.setSingleStranded(singleStranded);
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
		CompositionMap map = new CompositionMap(taxonMap,turnOffDestacking,turnOffDeDuping, useAllAlignments, filter);
		map.process();
		map.getNonStacked();
		map.setMapReader(mapReader);
		HashMap<String,ArrayList<Alignment>> list = map.getResultsMap();
		for(String key : list.keySet()){
			ArrayList<Alignment> als = list.get(key);
			String sequence = "";
			int count = 0;
			double length=0;
			double pIdent=0;
			double edit = 0;
			for(Alignment al : als) {
				numMatches++;	
				container.processAlignment(al);
				if(wantAlignments){
					alignments.add(al.getReadName());
					alignments.add(al.getText());
				}
				length += al.getReadLength();
				pIdent +=al.getPIdent();
				edit += al.getEditDistance();
				if( count == 0) {
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
				count++;
			}
			lengths.add((int)(length/=count));
			pIdents.add((pIdent/=count));
			distances.add((int)(edit/=count));
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
