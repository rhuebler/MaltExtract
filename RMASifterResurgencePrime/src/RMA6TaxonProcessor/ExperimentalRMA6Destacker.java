package RMA6TaxonProcessor;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.logging.Logger;

import NCBI_MapReader.NCBI_MapReader;
import RMAAlignment.Alignment;
import RMAAlignment.CompositionMap;
import behaviour.Filter;
import megan.data.IMatchBlock;
/**
 * RMA6 processor that automatically removes PCR duplicated and stacked reads Essentially it gets the alignment block 
 * from a read and than processes the best scoring percent of the alignment while sorting them by reference sequence
 * essentially first duplicate and stacking reads are marked and discarded and then authenticity criteria calculated
 * for the remaining reads
 * @author huebler
 * @params int id ,double pID, NCBI_MapReader reader, boolean v,Logger log, Logger warning, boolean reads,
 * double tp,int mL,boolean wantAls,boolean turnOffDestacking,boolean turnOffDeDuping,Filter behave
 */
public class ExperimentalRMA6Destacker extends RMA6TaxonProcessor {
	// initialize attributes
	protected boolean wantReads = false;
	protected boolean wantAlignments = false;
	protected boolean turnOffDestacking = false;
	protected boolean turnOffDeDuping = false;
	//construvters and set values
	public ExperimentalRMA6Destacker(Integer id, double pID, NCBI_MapReader reader, boolean v, Logger log, Logger warning,double tp,int mL, Filter behave) {
		super(id, pID, reader, v, log, warning,tp,mL, behave);
	}
	public ExperimentalRMA6Destacker(int id ,double pID, NCBI_MapReader reader,
			boolean v,Logger log, Logger warning, boolean reads,double tp,int mL,boolean wantAls,boolean turnOffDestacking,boolean turnOffDeDuping,Filter behave) {
		super(id,pID, reader, v, log, warning,tp,mL, behave);
		this.wantReads =reads;
		this.wantAlignments = wantAls;
		this.turnOffDestacking = turnOffDestacking;
		this.turnOffDeDuping = turnOffDeDuping;
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
		CompositionMap map = new CompositionMap(taxonMap,turnOffDestacking,turnOffDeDuping);
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
				String name = "";
				String readName = key;
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
	}//process
}
