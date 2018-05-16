package RMA6TaxonProcessor;
/**
 * The purpose of this class is to serve as mediator between the rma6File and the filters
 * it has different slots for each filters and retrieves the RMA6 match blocks from the file and than passes 
 * them into the filters to be processed and 
 * extract and keep some filters. This allows to run two or more filters in parallel.
 * @author huebler
 */

import java.util.logging.Level;
import java.util.logging.Logger;


import Analysis_16S_Data.RMA6_16S_AncientNodeProcessor;
import Analysis_16S_Data.RMA6_16S_NodeProcessor;
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
		private String fileName ="";
		private boolean downsample = true;
		//constructors
		public NodeProcessor(int id, double minPIdent, NCBI_MapReader reader, boolean v, Logger log, Logger warning,
				boolean reads, Filter behave, double minCompl,boolean alignment,boolean turnOffDestacking, 
				boolean turnOffDeDuping, boolean downsample) {
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
			this.downsample = downsample;
		}
		//getters
		
		public double getComplexity(String sequence){// calculate read complexity
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
		
		public int getTaxId(){
			return this.taxID;
		}
		public String getFileName(){
			return this.fileName;
		}
		//processing through node initialize taxon processors based on input parameters
		public void process(String inDir, String fileName, double topPercent, int maxLength){ 
			switch(behave) {
				case ANCIENT:
					ancientProcessor =  new ExperimentalRMA6AncientDestacker(taxID, minPIdent, mapReader, verbose, log, warning, wantReads, topPercent, maxLength,alignment,turnOffDestacking,turnOffDeDuping, behave);
					break;
				case NON:
					defaultProcessor = new ExperimentalRMA6Destacker(taxID, minPIdent, mapReader, verbose, log, warning, wantReads, topPercent, maxLength,alignment,turnOffDestacking,turnOffDeDuping, behave);
					break;
				case NON_ANCIENT:
					ancientProcessor = new ExperimentalRMA6AncientDestacker(taxID, minPIdent, mapReader, verbose, log, warning, wantReads, topPercent, maxLength,alignment,turnOffDestacking,turnOffDeDuping,  behave);
					defaultProcessor = new ExperimentalRMA6Destacker(taxID, minPIdent, mapReader, verbose, log, warning, wantReads, topPercent, maxLength,alignment,turnOffDestacking,turnOffDeDuping,  behave);
					break;
				case SRNA:
					ancientProcessor = new RMA6_16S_AncientNodeProcessor(taxID, minPIdent, mapReader, verbose, log, warning, wantReads, topPercent, maxLength,alignment,turnOffDestacking,turnOffDeDuping, behave);
					defaultProcessor = new RMA6_16S_NodeProcessor(taxID, minPIdent, mapReader, verbose, log, warning, wantReads, topPercent, maxLength,alignment,turnOffDestacking,turnOffDeDuping,  behave);
					break;
				default:
					System.err.println("Filter no longer supported");
					break;
					
			}
			this.taxName = getName(taxID);
			// use ReadsIterator to get all Reads assigned to MegantaxID and print top percent to file
			this.fileName =fileName;
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
				if(list.size()>10000 && downsample){
					warning.log(Level.WARNING,"For " + taxName + " in file "+fileName+ " downsampling was turned on");
				}
				IReadBlockIterator classIt  = new ReadBlockIterator(list, new ReadBlockGetterRMA6(rma6File, true, true, (float) 1.0,(float) 100.00,false,true));
				if(!classIt.hasNext()){ // check if reads are assigned to TaxID if not print to console and skip and set all values to default
					if(verbose){
						warning.log(Level.WARNING,"TaxID: " + taxName +  " not assigned in File " + fileName+"\n");
					}
				}else{
					if(verbose)
						log.log(Level.INFO,"Processing Taxon "+taxName+" in File " +fileName+"\n"
								+ " Number of Reads "+ list.size());
					int i = 0;
					while(classIt.hasNext()){// get Alignments and pass to filters
						IReadBlock current = classIt.next();
						if(current.getReadLength() <= maxLength || maxLength == 0){
							if(minComplexity<=getComplexity(current.getReadSequence())){
								IMatchBlock[] blocks = current.getMatchBlocks();
								switch(behave) {
								default:
									defaultProcessor.processMatchBlocks(blocks,current.getReadName(), current.getReadLength(), current.getReadSequence());
									ancientProcessor.processMatchBlocks(blocks, current.getReadName(), current.getReadLength(), current.getReadSequence());
									break;
								case ANCIENT:
									ancientProcessor.processMatchBlocks(blocks, current.getReadName(), current.getReadLength(), current.getReadSequence());
									break;
								case NON:
									defaultProcessor.processMatchBlocks(blocks,current.getReadName(), current.getReadLength(), current.getReadSequence());
									break;
								}
							}
						}// if  
						i++;
						if(i==10000 && downsample){
							break;
						}
					}// while
				classIt.close();
				rma6File.close();
				// process taxonprocessors
				switch(behave){
					default:
						ancientProcessor.process();
						defaultProcessor.process();
						ancientProcessor.clear();
						defaultProcessor.clear();
						break;
					case ANCIENT:	
						ancientProcessor.process();
						ancientProcessor.clear();
						break;
					case NON:
						defaultProcessor.process();
						defaultProcessor.clear();
						break;
						
				}
			}
		}catch(Exception e){
			warning.log(Level.SEVERE,mapReader.getNcbiIdToNameMap().get(taxID), e);	
		}
	}// void
}
