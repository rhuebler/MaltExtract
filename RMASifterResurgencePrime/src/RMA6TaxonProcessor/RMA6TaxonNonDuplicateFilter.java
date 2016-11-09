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
public class RMA6TaxonNonDuplicateFilter  extends RMA6TaxonProcessor{
	/**
	 * @param int ID, NCBI_MapReader reader, boolean verbose, Logger, log, Logger warning
	 * @return int numMatches, String readDistribution, HashMap EditDistance, HashMap Percent Identity
	 */ 
	public RMA6TaxonNonDuplicateFilter(int id ,double pID, NCBI_MapReader reader, boolean v,Logger log, Logger warning,double tp,int mL) {
		super(id,pID, reader, v, log, warning, tp, mL);
	}
	public void process(){
		CompositionMap map = new CompositionMap(taxonMap);
		map.process();
		map.markAllDuplicates();
		setReadDistribution(map);// first set ReadDistribution on Maximum ID
		taxonMap = map.getCompositionMap();
		for(int key : taxonMap.keySet()){
			for(Alignment entry : taxonMap.get(key)){
				if(!entry.isDuplicate()){
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
		setEditDistanceHistogram(distances);
		setPercentIdentityHistogram(pIdents);
	}	
	public void processMatchBlocks(IMatchBlock[] blocks, String readName, int readLength){ 
		IMatchBlock block = blocks[0];
		Alignment al = new Alignment();
		al.processText(block.getText().split("\n"));
		al.setReadName(readName);
		al.setPIdent(block.getPercentIdentity());
		if(minPIdent <= al.getPIdent()){
			if(!taxonMap.containsKey(block.getTaxonId())){
				ArrayList<Alignment> entry =new ArrayList<Alignment>();
				entry.add(al);
				taxonMap.put(block.getTaxonId(), entry);
			}else{
				ArrayList<Alignment> entry = taxonMap.get(block.getTaxonId());
				entry.add(al);
				taxonMap.put(block.getTaxonId(),entry);
				}
			}
		}//if 
	}
