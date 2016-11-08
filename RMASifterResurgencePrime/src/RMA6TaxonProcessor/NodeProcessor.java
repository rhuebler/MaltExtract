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
import strainMap.StrainMap;
import strainMap.StrainMisMatchContainer;
public class NodeProcessor{
		private RMA6TaxonDamageFilter ancientProcessor;
		private RMA6TaxonNonFilter defaultProcessor;
		private RMA6TaxonNonDuplicateFilter nonDuplicateProcessor;
		private RMA6TaxonAncientNonDuplicate ancientNonDuplicateProcessor;
		private String taxName;	
		private boolean wantReads = false;
		private NCBI_MapReader mapReader;
		private Integer taxID;
		private double minPIdent;
		private boolean verbose;
		private Logger log;
		private Logger warning;
		private Filter behave;
		
		public NodeProcessor(int id, double pID, NCBI_MapReader reader, boolean v, Logger log, Logger warning, Filter behave) {
			this.taxID = id;
			this.minPIdent = pID;
			this.mapReader = reader;
			this.verbose = v;
			this.log = log;
			this.warning = warning;
		}
		
		public NodeProcessor(int id, double pID, NCBI_MapReader reader, boolean v, Logger log, Logger warning,
				boolean reads, Filter behave) {
			this.taxID = id;
			this.minPIdent = pID;
			this.mapReader = reader;
			this.verbose = v;
			this.log = log;
			this.warning = warning;
			this.wantReads = reads;
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
				if(behave == Filter.NON_ANCIENT ||behave == Filter.ANCIENT ){
					
				}else if(behave == Filter.NON_ANCIENT ||behave == Filter.NON ){
					
				}else if(behave == Filter.ALL){
					
				}		
			}else{
if(behave == Filter.NON_ANCIENT ||behave == Filter.ANCIENT ){
					
				}else if(behave == Filter.NON_ANCIENT ||behave == Filter.NON ){
					
				}else if(behave == Filter.ALL){
					
				}else if(behave == Filter.NONDUPLICATES){
					
				}	
			}
			this.taxName = getName(taxID);
			int numMatches = 0;
			int ancientNumMatches = 0;
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
				if(!classIt.hasNext()){ // check if reads are assigned to TaxID if not print to console and skip
					if(verbose)
						warning.log(Level.WARNING,"TaxID: " + taxID +  " not assigned in File " + fileName+"\n");
					String s = taxName;
					for(int i = 0;i<=40;i++){
						s+="\t"+0;
					}
					defaultProcessor.setReadDistribution(new CompositionMap(new HashMap<Integer,ArrayList<Alignment>>()));
					defaultProcessor.setPercentIdentityHistogram(allPIdents);
					defaultProcessor.setEditDistanceHistogram(allDistances);
					defaultProcessor.setDamageLine(s);
					
					ancientProcessor.setDamageLine(s);
					ancientProcessor.setReadDistribution(new CompositionMap(new HashMap<Integer,ArrayList<Alignment>>()));
					ancientProcessor.setEditDistanceHistogram(ancientDistances);
					ancientProcessor.setPercentIdentityHistogram(ancientPIdents);
					
					nonDuplicateProcessor.setDamageLine(s);
					nonDuplicateProcessor.setReadDistribution(new CompositionMap(new HashMap<Integer,ArrayList<Alignment>>()));
					nonDuplicateProcessor.setEditDistanceHistogram(ancientDistances);
					nonDuplicateProcessor.setPercentIdentityHistogram(ancientPIdents);
					
					ancientNonDuplicateProcessor.setDamageLine(s);
					ancientNonDuplicateProcessor.setReadDistribution(new CompositionMap(new HashMap<Integer,ArrayList<Alignment>>()));
					ancientNonDuplicateProcessor.setEditDistanceHistogram(ancientDistances);
					ancientNonDuplicateProcessor.setPercentIdentityHistogram(ancientPIdents);
			}else{
				if(verbose)
					log.log(Level.INFO,"Processing Taxon "+taxName+" in File " +fileName); 
				HashMap<Integer, ArrayList<Alignment>> taxonMap = new HashMap<Integer,ArrayList<Alignment>>();
				HashMap<Integer, ArrayList<Alignment>> ancientMap = new HashMap<Integer,ArrayList<Alignment>>();
				int numReads = 0;
				int ancientNumReads = 0;
				while(classIt.hasNext()){
					IReadBlock current = classIt.next();
					boolean higher = false;
					int damage = 0;
					if(current.getReadLength() <= maxLength || maxLength == 0){
						IMatchBlock[] blocks=current.getMatchBlocks();
						int k=0;
						float topScore = current.getMatchBlock(0).getBitScore();
						double pIdent = 0;
						int editDistance=0;
						double ancientPIdent = 0;
						int ancientEditDistance=0;
						for(int i = 0; i< blocks.length;i++){
							if(blocks[i].getBitScore()/topScore < 1-topPercent){
								break;}
							
							numMatches++;
							Alignment al = new Alignment();
							al.processText(blocks[i].getText().split("\n"));
							al.setPIdent(blocks[i].getPercentIdentity());
							al.setReadName(current.getReadName());
							al.setReadLength(current.getReadLength());
							al.setAcessionNumber(blocks[i].getRefSeqId());	
							if(minPIdent <= al.getPIdent()){ // check for minPercentIdentity
								if(al.getFivePrimeDamage()){
									ancientNumMatches++;
									higher = true;
									//get mismatches
									ancientContainer.processAlignment(al);
									if(!ancientMap.containsKey(blocks[i].getTaxonId())){
										ArrayList<Alignment> entry =new ArrayList<Alignment>();
										entry.add(al);
										ancientMap.put(blocks[i].getTaxonId(), entry);
									}else{
										ArrayList<Alignment> entry = ancientMap.get(blocks[i].getTaxonId());
										entry.add(al);
										ancientMap.put(blocks[i].getTaxonId(),entry);
									}
									ancientEditDistance += al.getEditInstance();
									ancientPIdent += al.getPIdent();
									damage++;
									if(wantReads){
										String name = getName(blocks[i].getTaxonId());
										ancientLines.add(al.getReadName()+"\t"+"Length:\t"+al.getReadLength()+"\t");
										ancientLines.add(name+"\t"+al.getAccessionNumber()+"\t"+"Start:\t"+al.getStart()+"\t"+"End:\t"+al.getEnd());
										ancientLines.add("Q:\t"+al.getQuery());
										ancientLines.add("A:\t"+al.getAlignment());
										ancientLines.add("R:\t"+al.getReference()+"\n");
									}
								}
								//get mismatches
								allContainer.processAlignment(al);
								higher = true;
								pIdent += al.getPIdent();
								editDistance += al.getEditInstance();
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
							if(behave == Filter.NONDUPLICATES ||behave == Filter.ALL)
								break;
						}
						if(behave==Filter.NON || behave==Filter.NON_ANCIENT){
							if(higher){
								numReads++;
								allDistances.add(editDistance/k);
								allPIdents.add(pIdent/k);
							}
						}
						if(behave==Filter.ANCIENT || behave==Filter.NON_ANCIENT){
							if(higher&&damage != 0){
								ancientNumReads++;
								ancientDistances.add(ancientEditDistance/damage);
								ancientPIdents.add(ancientPIdent/damage);
							}
						}			
					}// if  
				}// while
					classIt.close();
					//set all non-filter information 
					if(behave==Filter.NON || behave==Filter.NON_ANCIENT){
						//set defaults
					CompositionMap map = new CompositionMap(taxonMap);
					map.process();
					StrainMap strain = new StrainMap(taxName,allContainer,numMatches);
					RMA6TaxonProcessor defaultProcessor = new RMA6TaxonProcessor();
					defaultProcessor.setDamageLine(strain.getLine());
					defaultProcessor.setNumberOfReads(numReads);
					defaultProcessor.setReadDistribution(map);
					defaultProcessor.setEditDistanceHistogram(allDistances);
					defaultProcessor.setPercentIdentityHistogram(allPIdents);
					this.defaultProcessor = defaultProcessor;
					}
					if(behave==Filter.ANCIENT || behave==Filter.NON_ANCIENT){
						//set ancient information
						CompositionMap aMap = new CompositionMap(ancientMap);
						aMap.process();
						StrainMap aStrain = new StrainMap(taxName,ancientContainer,ancientNumMatches);
						RMA6TaxonProcessor ancientProcessor = new RMA6TaxonProcessor();
						ancientProcessor.setDamageLine(aStrain.getLine());
						ancientProcessor.setNumberOfReads(ancientNumReads);
						ancientProcessor.setReadDistribution(aMap);
						ancientProcessor.setEditDistanceHistogram(ancientDistances);
						ancientProcessor.setPercentIdentityHistogram(ancientPIdents);
						ancientProcessor.setReads(ancientLines);
						this.ancientProcessor = ancientProcessor;
					}
					if(behave == Filter.NONDUPLICATES){
						computeOutputNonDefault(taxonMap,taxID);
					}
					if(behave == Filter.ALL){
						computeAncientNonDuplicates(ancientMap,taxID);
					}
					rma6File.close();
			     }//else
				}catch(Exception e){
					warning.log(Level.SEVERE,mapReader.getNcbiIdToNameMap().get(taxID), e);	
				
				}
			}// void
		private void computeAncientNonDuplicates(HashMap<Integer, ArrayList<Alignment>> taxonMap, int taxID){
			ArrayList<Integer> distances = new ArrayList<Integer>();
			ArrayList<Double> pIdents = new ArrayList<Double>();
			ArrayList<String> lines = new ArrayList<String>();
			StrainMisMatchContainer container = new StrainMisMatchContainer();
			int numMatches = 0;
			lines.add(taxName);
			CompositionMap map = new CompositionMap(taxonMap);
			map.process();
			map.markAllDuplicates();
			ancientNonDuplicateProcessor.setReadDistribution(map);
			taxonMap = map.getCompositionMap();
			int numReads=0;
			for(int key : taxonMap.keySet()){
				for(Alignment entry : taxonMap.get(key)){
					if(!entry.isDuplicate()){
						String name = getName(key);
						lines.add(entry.getReadName()+"\t"+"Length:\t"+entry.getReadLength()+"\t");
						lines.add(name+"\t"+entry.getAccessionNumber()+"\t"+"Start:\t"+entry.getStart()+"\t"+"End:\t"+entry.getEnd());
						lines.add("Q:\t"+entry.getQuery());
						lines.add("A:\t"+entry.getAlignment());
						lines.add("R:\t"+entry.getReference()+"\n");
						//get mismatches
						numMatches++;
						container.processAlignment(entry);
						pIdents.add(entry.getPIdent());
						distances.add(entry.getEditInstance());
						numReads++;
					}
					
				}
				
			}
			StrainMap strain = new StrainMap(taxName,container,numMatches);
			ancientNonDuplicateProcessor.setReads(lines);
			ancientNonDuplicateProcessor.setDamageLine(strain.getLine());
			ancientNonDuplicateProcessor.setNumberOfReads(numReads);
			ancientNonDuplicateProcessor.setEditDistanceHistogram(distances);
			ancientNonDuplicateProcessor.setPercentIdentityHistogram(pIdents);
			
		}	
		private void computeOutputNonDefault(HashMap<Integer, ArrayList<Alignment>> taxonMap, int taxID){
			ArrayList<Integer> distances = new ArrayList<Integer>();
			ArrayList<Double> pIdents = new ArrayList<Double>();
			ArrayList<String> lines = new ArrayList<String>();
			StrainMisMatchContainer container = new StrainMisMatchContainer();
			int numMatches = 0;
			lines.add(taxName);
			CompositionMap map = new CompositionMap(taxonMap);
			map.process();
			map.markAllDuplicates();
			nonDuplicateProcessor.setReadDistribution(map);
			taxonMap = map.getCompositionMap();
			int numReads=0;
			for(int key : taxonMap.keySet()){
				for(Alignment entry : taxonMap.get(key)){
					if(!entry.isDuplicate()){
						String name = getName(key);
						lines.add(entry.getReadName()+"\t"+"Length:\t"+entry.getReadLength()+"\t");
						lines.add(name+"\t"+entry.getAccessionNumber()+"\t"+"Start:\t"+entry.getStart()+"\t"+"End:\t"+entry.getEnd());
						lines.add("Q:\t"+entry.getQuery());
						lines.add("A:\t"+entry.getAlignment());
						lines.add("R:\t"+entry.getReference()+"\n");
						//get mismatches
						numMatches++;
						container.processAlignment(entry);
						pIdents.add(entry.getPIdent());
						distances.add(entry.getEditInstance());
						numReads++;
					}
					
				}
				
			}
			StrainMap strain = new StrainMap(taxName,container,numMatches);
			nonDuplicateProcessor.setReads(lines);
			nonDuplicateProcessor.setDamageLine(strain.getLine());
			nonDuplicateProcessor.setNumberOfReads(numReads);
			nonDuplicateProcessor.setEditDistanceHistogram(distances);
			nonDuplicateProcessor.setPercentIdentityHistogram(pIdents);
			
		}	
	}

