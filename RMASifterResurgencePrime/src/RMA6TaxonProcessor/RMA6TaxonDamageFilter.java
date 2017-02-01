package RMA6TaxonProcessor;

import java.util.ArrayList;
import java.util.logging.Logger;
import NCBI_MapReader.NCBI_MapReader;
import RMAAlignment.Alignment;
import RMAAlignment.CompositionMap;
import megan.data.IMatchBlock;
import strainMap.StrainMap;

/**
 * Extract all information from one Taxon and save the information in the specified slots to be retrieved
 * in RMA6Processor this taxon processor only processes reads that have at least one match that hints at c-> end substitutions 
 * it will be used to replace the non filtering taxon processor within RMA6Processor when filtering for damaged Reads 
 * @author huebler
 *
 */
public class RMA6TaxonDamageFilter extends RMA6TaxonProcessor{
	protected boolean wantReads = false;
	protected boolean wantAlignments = false;
	
	/**
	 * @param int ID, NCBI_MapReader reader, boolean verbose, Logger, log, Logger warning
	 * @return int numMatches, String readDistribution, HashMap EditDistance, HashMap Percent Identity
	 */ 
public RMA6TaxonDamageFilter(int id ,double pID, NCBI_MapReader reader,
		boolean v,Logger log, Logger warning, boolean reads,double tp,int mL,boolean wantAls) {
	super(id,pID, reader, v, log, warning, tp, mL);
	this.wantReads =reads;
	this.wantAlignments = wantAls;
}
public void processMatchBlocks(IMatchBlock[] blocks, String readName, int readLength,String sequence){ 
	int k=0;
	float topScore = blocks[0].getBitScore();
	double pIdent = 0;
	int editDistance = 0;
	int damage=0;
	boolean higher = false;
	for(int i = 0; i< blocks.length;i++){
		if(blocks[i].getBitScore()/topScore <= 1 - topPercent){
			break;}
		Alignment al = new Alignment();
		al.setSequence(sequence);
		al.setText(blocks[i].getText());
		al.processText();
		al.setPIdent(blocks[i].getPercentIdentity());
		al.setReadName(readName);
		al.setReadLength(readLength);
		al.setAcessionNumber(blocks[i].getRefSeqId());
		if(al.getFivePrimeDamage() && minPIdent <= al.getPIdent()){
			numMatches++;
			higher = true;
			//get mismatches
			container.processAlignment(al);
			if(!taxonMap.containsKey(blocks[i].getTaxonId())){
				ArrayList<Alignment> entry =new ArrayList<Alignment>();
				entry.add(al);
				taxonMap.put(blocks[i].getTaxonId(), entry);
			}else{
				ArrayList<Alignment> entry = taxonMap.get(blocks[i].getTaxonId());
				entry.add(al);
				taxonMap.put(blocks[i].getTaxonId(),entry);
			}
			editDistance += al.getEditDistance();
			pIdent += al.getPIdent();
			damage++;
			k++;
			if(wantAlignments){
				alignments.add(al.getText());
			}
		}
					
				}
				if(damage !=0 && higher){
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
					lengths.add(readLength);
					numOfReads++;
					distances.add(editDistance/k);
					pIdents.add(pIdent/k);
				}	
	}
		
	public void process(){
		StrainMap strain = new StrainMap(taxName,container,numMatches);
		CompositionMap map = new CompositionMap(taxonMap);
		map.process();
		setDamageLine(strain.getLine());
		setReadDistribution(map);
		setNumberOfReads(numOfReads);
		setEditDistanceHistogram(distances);
		setPercentIdentityHistogram(pIdents);
		setReads(lines);
		setAlignments(alignments);
		calculateReadLengthDistribution();
	 }
}// class 
