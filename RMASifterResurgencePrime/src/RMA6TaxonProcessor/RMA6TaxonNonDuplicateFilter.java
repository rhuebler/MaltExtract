package RMA6TaxonProcessor;

import java.util.ArrayList;
import java.util.logging.Logger;

import NCBI_MapReader.NCBI_MapReader;
import RMAAlignment.Alignment;
import RMAAlignment.CompositionMap;
import behaviour.Filter;
import megan.data.IMatchBlock;
import strainMap.StrainMap;
/**
 * NonDuplicate Filter for RMA6 Files
 * Due to technical constraints at the moment we only use the highest scoring BLast it when we remove duplicates at the moment
 * Runtime should be longer for duplicate removal due to passing all data one additional time 
 * @author huebler
 * @deprecated
 */
public class RMA6TaxonNonDuplicateFilter  extends RMA6TaxonProcessor{
	/**
	 * @param int ID, NCBI_MapReader reader, boolean verbose, Logger, log, Logger warning
	 * @return int numMatches, String readDistribution, HashMap EditDistance, HashMap Percent Identity
	 */ 
	protected boolean wantReads = false;
	protected boolean wantAlignment = false;
	public RMA6TaxonNonDuplicateFilter(int id ,double pID, NCBI_MapReader reader,
			boolean v,Logger log, Logger warning, boolean reads,double tp,int mL,boolean wantAlignments, Filter behave) {
		super(id,pID, reader, v, log, warning, tp, mL, behave);
		this.wantReads = reads;
		this.wantAlignment = wantAlignments;
	}
	public void process(){
		CompositionMap map = new CompositionMap(taxonMap);
		map.process();
		map.markAllDuplicates();
		processCompositionMap(map);// first set ReadDistribution on Maximum ID
		taxonMap = map.getCompositionMap();
		for(int key : taxonMap.keySet()){
			for(Alignment entry : taxonMap.get(key)){
				if(!entry.isDuplicate()){
					lengths.add(entry.getReadLength());
					//get mismatches
					if(wantAlignment){
						alignments.add(entry.getText());
					}
					if(wantReads){
						String name = "";
						if (!entry.getReadName().startsWith(">"))
	                        name = ">"+entry.getReadName();
						else
							name = entry.getReadName();
						lines.add(name);
	                    if (!name.endsWith("\n"))
	                        name += "\n";
	                    String readData = entry.getSequence();
	                    if (readData != null) {
	                        if (!readData.endsWith("\n"))
	                        	readData+=("\n");
	                    lines.add(readData);    
	                    }
					}
					numMatches++;
					container.processAlignment(entry);
					pIdents.add(entry.getPIdent());
					distances.add(entry.getEditDistance());
					numOfReads++;
				}
				
			}
			
		}
		StrainMap strain = new StrainMap(taxName,container,numMatches);
		setDamageLine(strain.getLine());
		setNumberOfReads(numOfReads);
		processCompositionMap(map);
		setEditDistanceHistogram(distances);
		setPercentIdentityHistogram(pIdents);
		calculateReadLengthDistribution();
		setAlignments(alignments);
		setReads(lines);
	}	
	public void processMatchBlocks(IMatchBlock[] blocks, String readName, int readLength,String sequence){ 
		IMatchBlock block = blocks[0];
		Alignment al = new Alignment();
		al.setSequence(sequence);
		al.setText(block.getText());
		al.processText();
		al.setReadName(readName);
		al.setPIdent(block.getPercentIdentity());
		al.setReadLength(readLength);
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
