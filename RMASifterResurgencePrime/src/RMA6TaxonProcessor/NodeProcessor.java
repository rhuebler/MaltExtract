package RMA6TaxonProcessor;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.logging.Level;
import java.util.logging.Logger;

import NCBI_MapReader.NCBI_MapReader;
import RMAAlignment.Alignment;
import RMAAlignment.CompositionMap;
import behaviour.Filter;
import jloda.util.ListOfLongs;
import megan.data.IMatchBlock;
import megan.data.IReadBlock;
import megan.data.IReadBlockIterator;
import megan.data.ReadBlockIterator;
import megan.rma6.ClassificationBlockRMA6;
import megan.rma6.RMA6File;
import megan.rma6.ReadBlockGetterRMA6;
public class NodeProcessor{
		private RMA6TaxonDamageFilter ancientProcessor;
		private RMA6TaxonNonFilter defaultProcessor;
		private RMA6TaxonNonDuplicateFilter nonDuplicateProcessor;
		private RMA6TaxonAncientNonDuplicate ancientNonDuplicateProcessor;
		private String taxName;	
		private boolean wantReads = false;
		private NCBI_MapReader mapReader;
		private Integer taxID;
		private boolean verbose;
		private Logger log;
		private Logger warning;
		private Filter behave;
		private double minPIdent;
		public NodeProcessor(int id,double minPIdent, NCBI_MapReader reader, boolean v, Logger log, Logger warning, Filter behave) {
			this.taxID = id;
			this.minPIdent = minPIdent;
			this.mapReader = reader;
			this.verbose = v;
			this.log = log;
			this.warning = warning;
			this.behave = behave;
		}
		
		public NodeProcessor(int id, double minPIdent, NCBI_MapReader reader, boolean v, Logger log, Logger warning,
				boolean reads, Filter behave) {
			this.taxID = id;
			this.minPIdent = minPIdent;
			this.mapReader = reader;
			this.verbose = v;
			this.log = log;
			this.warning = warning;
			this.wantReads = reads;
			this.behave = behave;
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
		public void process(String inDir, String fileName, double topPercent, int maxLength){ 
			if(wantReads){
				if(behave == Filter.ANCIENT){
					ancientProcessor = new RMA6TaxonDamageFilter(taxID, minPIdent, mapReader, verbose, log, log, wantReads, topPercent, maxLength);
				}else if(behave == Filter.NON){
					defaultProcessor = new RMA6TaxonNonFilter(taxID, minPIdent, mapReader, verbose, log, log, wantReads, topPercent, maxLength);
				}else if(behave == Filter.ALL){
					ancientNonDuplicateProcessor = new RMA6TaxonAncientNonDuplicate(taxID, minPIdent, mapReader, verbose, log, log, wantReads, topPercent, maxLength);
				}else if(behave == Filter.NON_ANCIENT){
					ancientProcessor = new RMA6TaxonDamageFilter(taxID, minPIdent, mapReader, verbose, log, log, wantReads, topPercent, maxLength);
					defaultProcessor = new RMA6TaxonNonFilter(taxID, minPIdent, mapReader, verbose, log, log, wantReads, topPercent, maxLength);
				}		
			}else{
				if(behave == Filter.ANCIENT ){
					ancientProcessor = new RMA6TaxonDamageFilter(taxID, minPIdent, mapReader, verbose, log, log, topPercent, maxLength);
				} 
				if(behave == Filter.NON ){
					defaultProcessor = new RMA6TaxonNonFilter(taxID, minPIdent, mapReader, verbose, log, log, topPercent, maxLength);
				}else if(behave == Filter.ALL){
					ancientNonDuplicateProcessor = new RMA6TaxonAncientNonDuplicate(taxID, minPIdent, mapReader, verbose, log, log, topPercent, maxLength);
				}else if(behave == Filter.NONDUPLICATES){
					nonDuplicateProcessor = new RMA6TaxonNonDuplicateFilter(taxID, minPIdent, mapReader, verbose, log, log, topPercent, maxLength);
				}else if(behave == Filter.NON_ANCIENT){
					ancientProcessor = new RMA6TaxonDamageFilter(taxID, minPIdent, mapReader, verbose, log, log, topPercent, maxLength);
					defaultProcessor = new RMA6TaxonNonFilter(taxID, minPIdent, mapReader, verbose, log, log, topPercent, maxLength);
				}	
			}
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
				IReadBlockIterator classIt  = new ReadBlockIterator(list, new ReadBlockGetterRMA6(rma6File, true, true, (float) 1.0,(float) 100.00,false,true));
				if(!classIt.hasNext()){ // check if reads are assigned to TaxID if not print to console and skip and set all values to default
					if(verbose)
						warning.log(Level.WARNING,"TaxID: " + taxID +  " not assigned in File " + fileName+"\n");
					String s = taxName;
					for(int i = 0;i<=40;i++){
						s+="\t"+0;
					}
					ArrayList<Double> pIdents = new ArrayList<Double>();
					ArrayList<Integer> distances = new ArrayList<Integer>();
					defaultProcessor.setReadDistribution(new CompositionMap(new HashMap<Integer,ArrayList<Alignment>>()));
					defaultProcessor.setPercentIdentityHistogram(pIdents);
					defaultProcessor.setEditDistanceHistogram(distances);
					defaultProcessor.setDamageLine(s);
					
					ancientProcessor.setDamageLine(s);
					ancientProcessor.setReadDistribution(new CompositionMap(new HashMap<Integer,ArrayList<Alignment>>()));
					ancientProcessor.setEditDistanceHistogram(distances);
					ancientProcessor.setPercentIdentityHistogram(pIdents);
					
					nonDuplicateProcessor.setDamageLine(s);
					nonDuplicateProcessor.setReadDistribution(new CompositionMap(new HashMap<Integer,ArrayList<Alignment>>()));
					nonDuplicateProcessor.setEditDistanceHistogram(distances);
					nonDuplicateProcessor.setPercentIdentityHistogram(pIdents);
					
					ancientNonDuplicateProcessor.setDamageLine(s);
					ancientNonDuplicateProcessor.setReadDistribution(new CompositionMap(new HashMap<Integer,ArrayList<Alignment>>()));
					ancientNonDuplicateProcessor.setEditDistanceHistogram(distances);
					ancientNonDuplicateProcessor.setPercentIdentityHistogram(pIdents);
			}else{
				if(verbose)
					log.log(Level.INFO,"Processing Taxon "+taxName+" in File " +fileName); 
				while(classIt.hasNext()){
					IReadBlock current = classIt.next();
					if(current.getReadLength() <= maxLength || maxLength == 0){
						IMatchBlock[] blocks=current.getMatchBlocks();
						if(behave == Filter.NON_ANCIENT ||behave == Filter.ANCIENT ){
							ancientProcessor.processMatchBlocks(blocks, current.getReadName(), current.getReadLength());
						} 
						if(behave == Filter.NON_ANCIENT ||behave == Filter.NON ){
							
							defaultProcessor.processMatchBlocks(blocks, current.getReadName(), current.getReadLength());
						}else if(behave == Filter.ALL){
							ancientNonDuplicateProcessor.processMatchBlocks(blocks, current.getReadName(), current.getReadLength());
						}else if(behave == Filter.NONDUPLICATES){
							nonDuplicateProcessor.processMatchBlocks(blocks, current.getReadName(), current.getReadLength());
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

