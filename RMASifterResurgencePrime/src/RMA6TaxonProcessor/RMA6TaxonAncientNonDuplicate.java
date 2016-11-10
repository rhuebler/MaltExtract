package RMA6TaxonProcessor;


import java.util.ArrayList;
import java.util.logging.Logger;
import NCBI_MapReader.NCBI_MapReader;
import RMAAlignment.Alignment;
import RMAAlignment.CompositionMap;
import megan.data.IMatchBlock;
import strainMap.StrainMap;
/**
 * NonDuplicate Filter for RMA6 Files
 * Due to technical constraints at the moment we only use the highest scoring BLast it when we remove duplicates at the moment
 * Runtime should be longer for duplicate removal due to passing all data one additional time 
 * @author huebler
 *
 */
public class RMA6TaxonAncientNonDuplicate  extends RMA6TaxonDamageFilter{
	/**
	 * @param int ID, NCBI_MapReader reader, boolean verbose, Logger, log, Logger warning
	 * @return int numMatches, String readDistribution, HashMap EditDistance, HashMap Percent Identity
	 */ 
	public RMA6TaxonAncientNonDuplicate(int id ,double pID, NCBI_MapReader reader,
			boolean v,Logger log, Logger warning,double tp,int mL) {
		super(id,pID, reader, v, log, warning,tp,mL);
	}
	public RMA6TaxonAncientNonDuplicate(int id ,double pID, NCBI_MapReader reader,
			boolean v,Logger log, Logger warning, boolean reads,double tp,int mL) {
		super(id,pID, reader, v, log, warning,reads,tp,mL);
	}
	
	public void process(){
		lines.add(taxName);
		CompositionMap map = new CompositionMap(taxonMap);
		map.process();
		map.markAllDuplicates();
		setReadDistribution(map);
		taxonMap = map.getCompositionMap();
		for(int key : taxonMap.keySet()){
			for(Alignment entry : taxonMap.get(key)){
				if(!entry.isDuplicate()){
					lengths.add(entry.getReadLength());
					if(wantReads){
						String name = getName(key);
						lines.add(entry.getReadName()+"\t"+"Length:\t"+entry.getReadLength()+"\t");
						lines.add(name+"\t"+entry.getAccessionNumber()+"\t"+"Start:\t"+entry.getStart()+"\t"+"End:\t"+entry.getEnd());
						lines.add("Q:\t"+entry.getQuery());
						lines.add("A:\t"+entry.getAlignment());
						lines.add("R:\t"+entry.getReference()+"\n");
					}
					//get mismatches
					numMatches++;
					container.processAlignment(entry);
					pIdents.add(entry.getPIdent());
					distances.add(entry.getEditInstance());
					numOfReads++;
				}
				
			}
			
		}
		StrainMap strain = new StrainMap(taxName,container,numMatches);
		setDamageLine(strain.getLine());
		setNumberOfReads(numOfReads);
		setReadDistribution(map);
		setEditDistanceHistogram(distances);
		setPercentIdentityHistogram(pIdents);
		setReads(lines);
		calculateReadLengthDistribution();
		
	}	
	public void processMatchBlocks(IMatchBlock[] blocks, String readName, int readLength){ 
		IMatchBlock block = blocks[0];
		Alignment al = new Alignment();
		al.processText(block.getText().split("\n"));
		al.setPIdent(block.getPercentIdentity());
		al.setReadName(readName);
		al.setReadLength(readLength);
		al.setAcessionNumber(block.getRefSeqId());
		if(al.getFivePrimeDamage() && minPIdent <= al.getPIdent()){
			if(!taxonMap.containsKey(block.getTaxonId())){
				ArrayList<Alignment> entry =new ArrayList<Alignment>();
				entry.add(al);
				taxonMap.put(block.getTaxonId(), entry);
			}else{
				ArrayList<Alignment> entry = taxonMap.get(block.getTaxonId());
				entry.add(al);
				taxonMap.put(block.getTaxonId(),entry);
				}
		}//5' damage	
	}//if 
	
}
