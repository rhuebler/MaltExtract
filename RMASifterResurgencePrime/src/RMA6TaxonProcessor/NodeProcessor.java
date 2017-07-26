package RMA6TaxonProcessor;
/**
 * The purpose of this class is to serve as mediator between the rma6File and the filters
 * it has different slots for each filters and passes the match blocks into the filters to be processed 
 * and keep some statistics 
 * @author huebler
 */

import java.util.ArrayList;
import java.util.Collections;
import java.util.logging.Level;
import java.util.logging.Logger;

import NCBI_MapReader.NCBI_MapReader;
import behaviour.Filter;
import jloda.util.DNAComplexityMeasure;
import jloda.util.ListOfLongs;
import megan.data.IMatchBlock;
import megan.data.IReadBlock;
import megan.data.IReadBlockIterator;
import megan.data.ReadBlockIterator;
import megan.rma6.ClassificationBlockRMA6;
import megan.rma6.RMA6File;
import megan.rma6.ReadBlockGetterRMA6;
public class NodeProcessor{
		private RMA6TaxonProcessor ancientProcessor;
		private RMA6TaxonProcessor defaultProcessor;
		private RMA6TaxonProcessor nonDuplicateProcessor;
		private RMA6TaxonProcessor ancientNonDuplicateProcessor;
		private String taxName;	
		private boolean wantReads = false;
		private NCBI_MapReader mapReader;
		private Integer taxID;
		private boolean verbose;
		private Logger log;
		private Logger warning;
		private Filter behave;
		private double minPIdent;
		private double minComplexity;
		private boolean alignment = false;
		private boolean turnOffDestacking = false;
		private boolean turnOffDeDuping = false;
		//constructors
		public NodeProcessor(int id, double minPIdent, NCBI_MapReader reader, boolean v, Logger log, Logger warning,
				boolean reads, Filter behave, double minCompl,boolean alignment,boolean turnOffDestacking, boolean turnOffDeDuping) {
			this.taxID = id;
			this.minPIdent = minPIdent;
			this.mapReader = reader;
			this.verbose = v;
			this.log = log;
			this.warning = warning;
			this.wantReads = reads;
			this.behave = behave;
			this.minComplexity = minCompl;
			this.alignment = alignment;
			this.turnOffDestacking = turnOffDestacking;
			this.turnOffDeDuping = turnOffDeDuping;
		}
		//getters
		// calculate complexity of read from Megan code
		public double getComplexity(String sequence){
			return DNAComplexityMeasure.getMinimumDNAComplexityWoottenFederhen(sequence);
		}
		private String getName(int taxId){
			String name;
			if(mapReader.getNcbiIdToNameMap().get(taxId) != null)
				name = mapReader.getNcbiIdToNameMap().get(taxId).replace(' ', '_');
			else
				name = "unassignedName";
			return name;
		}
		public RMA6TaxonProcessor getAncient(){
			return this.ancientProcessor;
			
		}
		public RMA6TaxonProcessor getDefault(){
			return  this.defaultProcessor;
			
		}
		public RMA6TaxonProcessor getNonDuplicate(){
			return this.nonDuplicateProcessor;
			
		}
		public RMA6TaxonProcessor getAncientNonDuplicate(){
			return  this.ancientNonDuplicateProcessor;
			
		}
		//processing
		public void process(String inDir, String fileName, double topPercent, int maxLength){ 
				if(behave == Filter.ANCIENT){
					ancientProcessor =  new ExperimentalRMA6AncientDestacker(taxID, minPIdent, mapReader, verbose, log, log, wantReads, topPercent, maxLength,alignment,turnOffDestacking,turnOffDeDuping, behave);
				}else if(behave == Filter.NON){
					defaultProcessor = new ExperimentalRMA6Destacker(taxID, minPIdent, mapReader, verbose, log, log, wantReads, topPercent, maxLength,alignment,turnOffDestacking,turnOffDeDuping, behave);
				}else if(behave == Filter.ALL){
					System.out.println("Filter no longer supported");
				}else if(behave == Filter.NON_ANCIENT){
					ancientProcessor = new ExperimentalRMA6AncientDestacker(taxID, minPIdent, mapReader, verbose, log, log, wantReads, topPercent, maxLength,alignment,turnOffDestacking,turnOffDeDuping,  behave);
					defaultProcessor = new ExperimentalRMA6Destacker(taxID, minPIdent, mapReader, verbose, log, log, wantReads, topPercent, maxLength,alignment,turnOffDestacking,turnOffDeDuping,  behave);
				}else if(behave == Filter.NONDUPLICATES ){
					System.out.println("Filter no longer supported");
				}		//TODO depreciated
			this.taxName = getName(taxID);
			// use ReadsIterator to get all Reads assigned to MegantaxID and print top percent to file
			try(RMA6File rma6File = new RMA6File(inDir+fileName, "r")){
				ListOfLongs list = new ListOfLongs();
				Long location = rma6File.getFooterSectionRMA6().getStartClassification("Taxonomy");
				if (location != null) {
				   ClassificationBlockRMA6 classificationBlockRMA6 = new ClassificationBlockRMA6("Taxonomy");
				   classificationBlockRMA6.read(location, rma6File.getReader());
				   if (classificationBlockRMA6.getSum(taxID) > 0) {
					   classificationBlockRMA6.readLocations(location, rma6File.getReader(), taxID, list);
				   }
				 }
				// Downsample list if necessary and write log 
				if(list.size()>100000){
					warning.log(Level.WARNING,"For " + taxName + " in file "+fileName+ " downsampling was turned on");
					ArrayList<Long> longList = new ArrayList<Long>();
					for(int i = 0; i<list.size(); i++)
						longList.add(list.get(i));
					Collections.shuffle(longList);
					list.clear();
					for(long l : longList.subList(0, 100000))
						list.add(l);
				}
				IReadBlockIterator classIt  = new ReadBlockIterator(list, new ReadBlockGetterRMA6(rma6File, true, true, (float) 1.0,(float) 100.00,false,true));
				if(!classIt.hasNext()){ // check if reads are assigned to TaxID if not print to console and skip and set all values to default
					if(verbose)
						warning.log(Level.WARNING,"TaxID: " + taxID +  " not assigned in File " + fileName+"\n");
				}else{
					if(verbose)
						log.log(Level.INFO,"Processing Taxon "+taxName+" in File " +fileName); 
					while(classIt.hasNext()){
						IReadBlock current = classIt.next();
						if(current.getReadLength() <= maxLength || maxLength == 0){
							if(minComplexity<=getComplexity(current.getReadSequence())){
								IMatchBlock[] blocks=current.getMatchBlocks();
								if(behave == Filter.NON_ANCIENT ||behave == Filter.ANCIENT ){
									ancientProcessor.processMatchBlocks(blocks, current.getReadName(), current.getReadLength(),current.getReadSequence());
								} 
								if(behave == Filter.NON_ANCIENT ||behave == Filter.NON ){
									defaultProcessor.processMatchBlocks(blocks, current.getReadName(), current.getReadLength(), current.getReadSequence());
								}else if(behave == Filter.ALL){
									ancientNonDuplicateProcessor.processMatchBlocks(blocks, current.getReadName(), current.getReadLength(), current.getReadSequence());
								}else if(behave == Filter.NONDUPLICATES){
									nonDuplicateProcessor.processMatchBlocks(blocks, current.getReadName(), current.getReadLength(), current.getReadSequence());
								}
							}
						}// if  
					}// while
				classIt.close();
				rma6File.close();
				if(behave == Filter.ANCIENT ){
					ancientProcessor.process();
				}else if(behave == Filter.NON ){
					defaultProcessor.process();
				}else if(behave == Filter.ALL){
					ancientNonDuplicateProcessor.process();
				}else if(behave == Filter.NONDUPLICATES){
					nonDuplicateProcessor.process();
				}else if(behave == Filter.NON_ANCIENT) {
					ancientProcessor.process();
					defaultProcessor.process();
				}	
			}
		}catch(Exception e){
			warning.log(Level.SEVERE,mapReader.getNcbiIdToNameMap().get(taxID), e);	
		}
	}// void
		
}

