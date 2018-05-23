package RMA6TaxonProcessor;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.concurrent.ConcurrentLinkedDeque;
import java.util.logging.Logger;

import NCBI_MapReader.NCBI_MapReader;
import RMAAlignment.Alignment;
import RMAAlignment.CompositionMap;
import behaviour.Filter;
/**
 * RMA6 child of RMA6TaxonProcessor used for RMA6Crawler. Has an extra function that allows to go through
 * a list RMA6 Alignmet projects and retrieve statistics for them all
 * @params int id ,double pID, NCBI_MapReader reader, boolean v,Logger log, Logger warning, boolean reads, 
 * double tp,int mL,boolean wantAls,boolean turnOffDestacking,boolean turnOffDeDuping,Filter behave
 */
public class MatchProcessorCrawler extends RMA6TaxonProcessor {
	// initialize attributes
	protected boolean wantReads = false;
	protected boolean wantAlignments = false;
	protected boolean turnOffDestacking = false;
	protected boolean turnOffDeDuping = false;
	//constructors and set values

	public MatchProcessorCrawler(int id ,double pID, NCBI_MapReader reader,boolean v,Logger log, Logger warning, boolean reads,double tp,int mL,boolean wantAls, boolean turnOffDestacking,boolean turnOffDeDuping, Filter behave, boolean useAllAlignments) {
		//intialize Nodeanalyzer
		super(id ,pID, reader, v, log, warning, reads, tp, mL, wantAls, turnOffDestacking, turnOffDeDuping, behave, useAllAlignments);
	}
	//process each Matchblock
	public void processDLQlist(ConcurrentLinkedDeque<Alignment> concurrentLinkedDeque) {
		for(Alignment al: concurrentLinkedDeque) {
			processMatchBlock(al);
		}
	}
	public void processMatchBlock(Alignment al){
			//System.out.println("Here");
				if(minPIdent <= al.getPIdent()){ // check for minPercentIdentity
					if(wantReads){//add Read
						String sequence = "";
						String name = "";
						String readName =al.getReadName();
						sequence = al.getSequence();
						if (!readName.startsWith(">"))
							name = ">"+readName;
						else
							name = readName;
						lines.add(name);
						if (!name.endsWith("\n"))
							name += "\n";
						String readData = sequence;
						if (readData != null) {
							if (!readData.endsWith("\n"))
								readData+=("\n");
							lines.add(readData);    
						}
					}
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
	public void process(){ 
		//analyze collected data 
		CompositionMap map = new CompositionMap(taxonMap,turnOffDestacking,turnOffDeDuping,useAllAlignments,filter);
		map.setMapReader(mapReader);
		map.process();
		map.getNonStacked();
		
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
		setTurnedOn(map.wasTurnedOn());
		calculateReadLengthDistribution();
		map = null;
	}
}
