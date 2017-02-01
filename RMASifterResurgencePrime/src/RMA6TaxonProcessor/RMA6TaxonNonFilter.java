package RMA6TaxonProcessor;

import java.util.ArrayList;
import java.util.logging.Logger;

import NCBI_MapReader.NCBI_MapReader;
import RMAAlignment.Alignment;
import RMAAlignment.CompositionMap;
import megan.data.IMatchBlock;
import strainMap.StrainMap;
/**
 * Default filter that just outputs additional proofs without any filtering
 * @author huebler
 *
 */
public class RMA6TaxonNonFilter  extends RMA6TaxonProcessor{
	protected boolean wantReads = false;
	protected boolean wantAlignments = false;
	public RMA6TaxonNonFilter(Integer id, double pID, NCBI_MapReader reader, boolean v, Logger log, Logger warning,double tp,int mL) {
		super(id, pID, reader, v, log, warning,tp,mL);
	}
	public RMA6TaxonNonFilter(int id ,double pID, NCBI_MapReader reader,
			boolean v,Logger log, Logger warning, boolean reads,double tp,int mL,boolean wantAls) {
		super(id,pID, reader, v, log, warning,tp,mL);
		this.wantReads =reads;
		this.wantAlignments = wantAls;
	}
	
	public void processMatchBlocks(IMatchBlock[] blocks, String readName, int readLength, String sequence){
		lengths.add(readLength);
		int k=0;
		float topScore = blocks[0].getBitScore();
		double pIdent = 0;
		int editDistance = 0;
		boolean higher = false;
			for(int i = 0; i< blocks.length;i++){
				if(blocks[i].getBitScore()/topScore < 1-topPercent){
					break;}
				numMatches++;			
				Alignment al = new Alignment();
				al.setText(blocks[i].getText());
				al.processText();
				al.setPIdent(blocks[i].getPercentIdentity());
				al.setReadName(readName);
				al.setReadLength(readLength);
				al.setAcessionNumber(blocks[i].getRefSeqId());	
				if(minPIdent <= al.getPIdent()){ // check for minPercentIdentity
					if(wantAlignments){
						alignments.add(al.getText());
					}
								//get mismatches
				container.processAlignment(al);
				higher = true;
				pIdent += al.getPIdent();
				editDistance += al.getEditDistance();
				if(!taxonMap.containsKey(blocks[i].getTaxonId())){
					ArrayList<Alignment> entry =new ArrayList<Alignment>();
					entry.add(al);
					taxonMap.put(blocks[i].getTaxonId(), entry);
				}else{
					ArrayList<Alignment> entry = taxonMap.get(blocks[i].getTaxonId());
					entry.add(al);
					taxonMap.put(blocks[i].getTaxonId(),entry);
					}
				k++;
			}
		}
			if(higher){
				numOfReads++;
				distances.add(editDistance/k);
				pIdents.add(pIdent/k);
				if(wantReads){
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
}		
	public void process(){ 
		
		CompositionMap map = new CompositionMap(taxonMap);
		map.process();
		StrainMap strain = new StrainMap(taxName,container,numMatches);
		setDamageLine(strain.getLine());
		setNumberOfReads(numOfReads);
		setNumMatches(numMatches);
		setReadDistribution(map);
		setEditDistanceHistogram(distances);
		setPercentIdentityHistogram(pIdents);
		setReads(lines);
		setAlignments(alignments);
		calculateReadLengthDistribution();
	}//process
		
}
