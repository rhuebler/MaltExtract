package utility;

import java.io.IOException;
import java.util.Set;
import java.util.concurrent.ConcurrentHashMap;
import java.util.logging.Level;
import java.util.logging.Logger;

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
	private ConcurrentHashMap<String, MatchProcessorCrawler> concurrentMap;
	private String filePath="";
	int id =0;
	Logger warning;
	boolean wantReads;
	double topPercent;
	String speciesName;
	public NodeMatchSorter(ConcurrentHashMap<String, MatchProcessorCrawler> concurrentMap,String filePath, int id, Logger warning,boolean wantReads,
			String speciesName){
		this.concurrentMap = concurrentMap;
		this.filePath = filePath;
		this.id = id;
		this.warning = warning;
		this.speciesName = speciesName;
	}
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
					if(blocks[i].getBitScore()/topScore < 1-topPercent){
						break;}		
					Alignment al = new Alignment();
					al.setText(blocks[i].getText());
					al.processText();
					al.setTaxID(blocks[i].getTaxonId());
					al.setPIdent(blocks[i].getPercentIdentity());
					al.setReadName(readName);
					al.setReadLength(readLength);
					al.setAcessionNumber(blocks[i].getTextFirstWord());	
					al.setSequence(readSequence);
					
				}
			}	
			classIt.close();
			rma6File.close();
		}catch(IOException io){
			warning.log(Level.SEVERE,"Can not locate or read File" ,io);
		}
	}
}
