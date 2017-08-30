package RMA6TaxonProcessor;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.logging.Logger;

import NCBI_MapReader.NCBI_MapReader;
import RMAAlignment.Alignment;
import RMAAlignment.CompositionMap;
import behaviour.Filter;
import megan.data.IMatchBlock;
import strainMap.StrainMap;
/**
 * RMA6 processor that automatically removes PCR duplicated and stacked reads
 * @author huebler
 *
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
	//process each Marchblock
	public void processMatchBlocks(IMatchBlock[] blocks, String readName, int readLength, String sequence){
		originalNumberOfReads++;
		System.out.println(blocks.length);
		float topScore = blocks[0].getBitScore();
			for(int i = 0; i< blocks.length;i++){
				if(blocks[i].getBitScore()/topScore < 1-topPercent){
					break;}		
			
				Alignment al = new Alignment();
				al.setText(blocks[i].getText());
				al.processText();
				al.setPIdent(blocks[i].getPercentIdentity());
				System.out.println((int)blocks[i].getPercentIdentity());//TODO Remove
				System.out.println(blocks[i].getText());
				al.setReadName(readName);
				al.setReadLength(readLength);
				al.setAcessionNumber(blocks[i].getTextFirstWord().split("\\|")[0].substring(1));
			
				al.setSequence(sequence);
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
		HashMap<String,ArrayList<Alignment>> list =map.getResultsMap();
		for(String key : list.keySet()){
			int k = 0;
			double pIdent = 0;
			int editDistance = 0;
			String sequence = "";
			
			for(Alignment al : list.get(key)){
				if(wantAlignments&&k==0)
					
				numMatches++;	
				pIdent+= al.getPIdent();
				editDistance += al.getEditDistance();
				container.processAlignment(al);
				if(wantAlignments){
					alignments.add(al.getReadName());
					alignments.add(al.getText());}
				if(k==0){
					lengths.add(al.getReadLength());
					if(wantReads)
						sequence=al.getSequence();
				}	
				k++;
			}
			if(k!=0){
				pIdents.add(pIdent/k);
				distances.add(Math.round(editDistance/k));
				if(wantReads){
					String readName =key;
					String name = "";
					if (!readName.startsWith(">"))
	                    name = ">"+readName;
					else
						name = readName;
					lines.add(name+"\t"+Math.round(editDistance/k));
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
		}
		StrainMap strain = new StrainMap(taxName,container,numMatches);
		setOriginalNumberOfAlignments(originalNumberOfAlignments);
		setOriginalNumberOfReads(originalNumberOfReads);
		setDamageLine(strain.getLine());
		setNumMatches(numMatches);
		setNumberOfReads(list.keySet().size());
		processCompositionMap(map);
		setEditDistanceHistogram(distances);
		setPercentIdentityHistogram(pIdents);
		setReads(lines);
		setAlignments(alignments);
		setTurnedOn(map.wasTurnedOn());
		calculateReadLengthDistribution();
		strain = null;
		
	}//process
}
