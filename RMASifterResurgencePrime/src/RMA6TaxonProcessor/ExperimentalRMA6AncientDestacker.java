package RMA6TaxonProcessor;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.logging.Logger;

import NCBI_MapReader.NCBI_MapReader;
import RMAAlignment.Alignment;
import RMAAlignment.CompositionMap;
import megan.data.IMatchBlock;
import strainMap.StrainMap;

public class ExperimentalRMA6AncientDestacker extends RMA6TaxonProcessor {
	protected boolean wantReads = false;
	protected boolean wantAlignments = false;
	
	public ExperimentalRMA6AncientDestacker(Integer id, double pID, NCBI_MapReader reader, boolean v, Logger log, Logger warning,double tp,int mL) {
		super(id, pID, reader, v, log, warning,tp,mL);
	}
	public ExperimentalRMA6AncientDestacker(int id ,double pID, NCBI_MapReader reader,
			boolean v,Logger log, Logger warning, boolean reads,double tp,int mL,boolean wantAls) {
		super(id,pID, reader, v, log, warning,tp,mL);
		this.wantReads =reads;
		this.wantAlignments = wantAls;
	}
	
	public void processMatchBlocks(IMatchBlock[] blocks, String readName, int readLength, String sequence){
		originalNumberOfReads++;
		float topScore = blocks[0].getBitScore();
			for(int i = 0; i< blocks.length;i++){
				if(blocks[i].getBitScore()/topScore < 1-topPercent){
					break;}		
				Alignment al = new Alignment();
				al.setText(blocks[i].getText());
				al.processText();
				al.setPIdent(blocks[i].getPercentIdentity());
				al.setReadName(readName);
				al.setReadLength(readLength);
				al.setAcessionNumber(blocks[i].getRefSeqId());	
				
				if(al.getFivePrimeDamage()&&minPIdent <= al.getPIdent()){ // check for minPercentIdentity
					originalNumberOfAlignments++;
								//get mismatches
					if(!taxonMap.containsKey(blocks[i].getTaxonId())){
						ArrayList<Alignment> entry =new ArrayList<Alignment>();
						entry.add(al);
						taxonMap.put(blocks[i].getTaxonId(), entry);
					}else{
						ArrayList<Alignment> entry = taxonMap.get(blocks[i].getTaxonId());
						entry.add(al);
						taxonMap.put(blocks[i].getTaxonId(),entry);
					}
				}
			}
	}		
	public void process(){ 
		
		CompositionMap map = new CompositionMap(taxonMap);
		map.process();
		map.getNonStacked();
		HashMap<String,ArrayList<Alignment>> list =map.getResultsMap();
		for(String key : list.keySet()){
			int k = 1;
			double pIdent = 0;
			int editDistance = 0;
			String sequence = "";
			for(Alignment al : list.get(key)){
				numMatches++;	
				pIdent+= al.getPIdent();
				editDistance = al.getEditDistance();
				container.processAlignment(al);
				if(wantAlignments)
					alignments.add(al.getAlignment());
				if(k==1){
					lengths.add(al.getReadLength());
					if(wantReads)
						sequence=al.getSequence();
				}	
				k++;
			}
			pIdents.add(pIdent/k);
			distances.add(editDistance/k);
			if(wantReads){
				String readName =key;
				String name = "";
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
		
		StrainMap strain = new StrainMap(taxName,container,numMatches);
		setOriginalNumberOfAlignments(originalNumberOfAlignments);
		setOriginalNumberOfReads(originalNumberOfReads);
		setDamageLine(strain.getLine());
		setNumMatches(numMatches);
		setNumberOfReads(list.keySet().size());
		setReadDistribution(map);
		setEditDistanceHistogram(distances);
		setPercentIdentityHistogram(pIdents);
		setReads(lines);
		setAlignments(alignments);
		calculateReadLengthDistribution();
	}//process
}
