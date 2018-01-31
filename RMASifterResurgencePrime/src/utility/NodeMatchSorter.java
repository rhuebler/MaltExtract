package utility;
/**
 * @defunct
 */
import java.io.IOException;
import java.util.HashMap;
import java.util.logging.Level;
import java.util.logging.Logger;

import NCBI_MapReader.NCBI_MapReader;
import RMA6TaxonProcessor.MatchProcessorCrawler;
import RMAAlignment.Alignment;
import jloda.util.ListOfLongs;
import megan.data.IMatchBlock;
import megan.data.IReadBlock;
import megan.data.IReadBlockIterator;
import megan.data.ReadBlockIterator;
import megan.rma6.ClassificationBlockRMA6;
import megan.rma6.RMA6File;
import megan.rma6.ReadBlockGetterRMA6;

public class NodeMatchSorter {
	private HashMap<String, MatchProcessorCrawler> concurrentMap= new HashMap<String, MatchProcessorCrawler>();
	private String filePath="";
	int id =0;
	Logger log;
	Logger warning;
	boolean wantReads;
	double topPercent;
	String speciesName;
	NCBI_MapReader mapReader;
	//Constructor
	public HashMap<String, MatchProcessorCrawler> returnCHashMap(){
		return concurrentMap;
	}
	public NodeMatchSorter(String filePath, int id,Logger log, Logger warning,boolean wantReads,
			String speciesName,NCBI_MapReader mapReader){
		this.filePath = filePath;
		this.id = id;
		this.warning = warning;
		this.speciesName = speciesName;
		this.mapReader=mapReader;
		this.log = log;
		this.wantReads =wantReads;
	}
	private String getName(int taxId){
		String name;
		if(mapReader.getNcbiIdToNameMap().get(taxId) != null)
			name = mapReader.getNcbiIdToNameMap().get(taxId);
		else
			name = "unassignedName";
		return name;
	}
	//process a whole Node to resort the matches to different strains
	public void processNode(){
		try(RMA6File rma6File = new RMA6File(filePath, "r")){
			ListOfLongs list = new ListOfLongs();
			Long location = rma6File.getFooterSectionRMA6().getStartClassification("Taxonomy");
			if (location != null) {
			   ClassificationBlockRMA6 classificationBlockRMA6 = new ClassificationBlockRMA6("Taxonomy");
			   classificationBlockRMA6.read(location, rma6File.getReader());
			   if (classificationBlockRMA6.getSum(id) > 0) {
				   classificationBlockRMA6.readLocations(location, rma6File.getReader(), id, list);
			   }
			 }
			
			IReadBlockIterator classIt  = new ReadBlockIterator(list, new ReadBlockGetterRMA6(rma6File, true, true, (float) 1.0,(float) 100.00,false,true));
			
			// iterate through all nodes and store information in strain Map to process and retrieve at later use
			while(classIt.hasNext()){
				IReadBlock current = classIt.next();
				IMatchBlock[] blocks = current.getMatchBlocks();
				String readName = current.getReadName();
				String readSequence = current.getReadSequence();
				int readLength = current.getReadLength();
				float topScore = blocks[0].getBitScore();
				for(int i = 0; i< blocks.length;i++){
					if((blocks[i].getBitScore()/topScore) < 1-topPercent){
						break;}	
					String name = getName(blocks[i].getTaxonId());
					if(name.contains(speciesName)){
						Alignment al = new Alignment();
						al.setText(blocks[i].getText());
						al.processText();
						al.setTaxID(blocks[i].getTaxonId());
						al.setPIdent(blocks[i].getPercentIdentity());
						al.setReadName(readName);
						al.setReadLength(readLength);
						al.setAcessionNumber(blocks[i].getTextFirstWord());	
						al.setSequence(readSequence);
						al.setTopAlignment(true);
						if(concurrentMap.containsKey(name)){
							MatchProcessorCrawler mpc = concurrentMap.get(name);
							mpc.processMatchBlock(al);
							concurrentMap.replace(name, mpc);
						}
						else{
							MatchProcessorCrawler mpc = new MatchProcessorCrawler(id,topPercent,mapReader,false,log,warning,wantReads,0.01,0,false,true,false,behaviour.Filter.CRAWL);
							mpc.processMatchBlock(al);
							concurrentMap.put(name, mpc);
						}
						break;
					}
				}
			}	
			classIt.close();
			rma6File.close();
		}catch(IOException io){
			warning.log(Level.SEVERE,"Can not locate or read File" ,io);
		}
	}
}
