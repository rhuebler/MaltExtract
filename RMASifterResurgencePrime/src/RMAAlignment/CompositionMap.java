package RMAAlignment;
import java.util.ArrayList;
import java.util.HashMap;
import NCBI_MapReader.NCBI_MapReader;
import behaviour.Filter;
/**
 * Composition Map is used to process a Hashmap that maps the ALignment objects to their 
 * specific references so that statistics concerning read alignments can be calculated correctly
 * It can also mark PCR duplicates
 * @author huebler
 */
// set Composition Map
public class CompositionMap {
	private boolean useAllAlignments = false;
	private boolean turnOffDestacking= false;
	private boolean turnOffDeDupping = false;
	private Filter behave = Filter.NON_ANCIENT;
//	public CompositionMap(HashMap<Integer,ArrayList<Alignment>> map){
//		setCompositionMap(map);
//	}
	public CompositionMap(HashMap<Integer, HashMap<String, ArrayList<Alignment>>> taxonMap,boolean turnOffDestacking, boolean turnOffDeDupping, boolean useAllAlignments){
		setCompositionMap(taxonMap);
		this.turnOffDestacking = turnOffDestacking;
		this.turnOffDeDupping = turnOffDeDupping;
		this.useAllAlignments = useAllAlignments;
	}
	public CompositionMap(HashMap<Integer, HashMap<String, ArrayList<Alignment>>> taxonMap,boolean turnOffDestacking, boolean turnOffDeDupping, boolean useAllAlignments,Filter behave){
		setCompositionMap(taxonMap);
		this.turnOffDestacking = turnOffDestacking;
		this.turnOffDeDupping = turnOffDeDupping;
		this.behave=behave;
		this.useAllAlignments = useAllAlignments;
	}
	// initialize values
private HashMap<Integer, HashMap<String, ArrayList<Alignment>>> compositionMap;// hashMap of ReferenceID to List of start positions
private HashMap<String,ArrayList<Alignment>> resultsMap = new HashMap<String,ArrayList<Alignment>>();// hashMap of ReferenceID to List of start positions
private int maxID=0;
private String maxReference;
private ArrayList<Double> generalStatistics;
private HashMap<Integer,Integer> coverageHistogram;
private ArrayList<String> coveragePositions;
private ArrayList<NOAOR> stackedSizes = new ArrayList<NOAOR>();
private int refLength = 0;
private boolean turnedOn = true;
private NCBI_MapReader mapReader;
private int maxSize;
//getter
public ArrayList<String> getCoveragePositions(){
	return this.coveragePositions;
}
public boolean wasTurnedOn(){
	return this.turnedOn;
}
public HashMap<String,ArrayList<Alignment>> getResultsMap(){
	return this.resultsMap;
}
public int getReferenceLength(){
	return this.refLength;
}
public int getMaxID(){
	return this.maxID;}

public HashMap<Integer, HashMap<String, ArrayList<Alignment>>> getCompositionMap(){
	return this.compositionMap;}
public ArrayList<Double> getGenaralStatistics(){
	return this.generalStatistics;
}
public HashMap<Integer,Integer> getConverageHistogram(){
	return this.coverageHistogram;
}
//setters
public void setMapReader(NCBI_MapReader mapReader){
	this.mapReader = mapReader;
}
private void setCompositionMap(HashMap<Integer, HashMap<String, ArrayList<Alignment>>> taxonMap){// technically not required
	this.compositionMap = taxonMap;}

private void setMaxID(int i){
	this.maxID = i;
}
private void setMaxReference(String reference){
	this.maxReference = reference;
}
//utility
private HashMap<String,ArrayList<Alignment>> markDuplicates(HashMap<String,ArrayList<Alignment>> map){ // that seems to be correct may be the ordering works better and therefore more stacking reads can be found
	AlignmentComparator comp = new AlignmentComparator(); // should theoretically sort my alignments according to start positions
	for(String key : map.keySet()){
		ArrayList<Alignment> inList = map.get(key);
		inList.sort(comp);
		int i =0;
		HashMap<String, Integer> timesViewed =  new HashMap<String, Integer>();
		while(i < inList.size() - 1){// array size is 47 last 46
			
			Alignment current = inList.get(i);
			if(timesViewed.containsKey(current.getQuery())){
				int k=timesViewed.get(current.getQuery());
				timesViewed.replace(current.getQuery(),k );
			}else{
				timesViewed.put(current.getQuery(), 1);
			}
			
			Alignment next = inList.get(i+1);
			int cStart = 0;
			int cEnd = 0;
			int nStart = 0;
			int nEnd = 0;
			if(!current.isReversed() == !next.isReversed()){
				cStart = current.getStart();
				cEnd = current.getEnd();
				nStart = next.getStart();
				nEnd = next.getEnd();
			}else if(current.isReversed() == !next.isReversed()){
				cStart = current.getEnd();
				cEnd = current.getStart();
				nStart = next.getStart();
				nEnd = next.getEnd();
			}else if(!current.isReversed() == next.isReversed()){
				cStart = current.getStart();
				cEnd = current.getEnd();
				nStart = next.getEnd();
				nEnd = next.getStart();
			}else{
				cStart = current.getEnd();
				cEnd = current.getStart();
				nStart = next.getEnd();
				nEnd = next.getStart();
			}
			if( cStart == nStart && cEnd == nEnd && current.getReferenceLength() == next.getReferenceLength()){
				if(timesViewed.containsKey(next.getQuery())||current.getQuery().equals(next.getQuery())){
				inList.get(i+1).setDuplicate(true);
				}
			}
				i++;
		}
		map.replace(key, inList);
	}
	return map;
}
public void markReferenceDuplicates(){
	this.compositionMap.put(getMaxID(), markDuplicates(this.compositionMap.get(getMaxID())));
}
public void markAllDuplicates(){
	for(int key : this.compositionMap.keySet()){
	this.compositionMap.put(key, markDuplicates(this.compositionMap.get(key)));
	}
}
public void calculateStatistics(){
	AlignmentStatistics stats = new AlignmentStatistics(compositionMap.get(getMaxID()).get(maxReference),turnOffDestacking, turnOffDeDupping, behave);
	stats.calculateStatistics();
	this.coverageHistogram = stats.getConverageHistogram();
	this.generalStatistics = stats.getGenaralStatistics();
	this.refLength = stats.getLength();
	this.coveragePositions = stats.getCoveragePositions();
	this.maxSize = stats.getDestackedList().size();
}

// process composition and find taxon with maximum number of start positions
public void getNonStacked(){
	//get nonstacked Reads
	HashMap<String,ArrayList<Alignment>> map = new HashMap<String,ArrayList<Alignment>>();
	for(int key : compositionMap.keySet()){
		HashMap<String,Integer> nonStackedOnReference = new HashMap<String,Integer>();
		HashMap<String,ArrayList<Alignment>> rMap = compositionMap.get(key);
		for(String reference : rMap.keySet()){
			if(turnOffDestacking){
				turnedOn = false;
				ArrayList<Alignment>list=rMap.get(reference);
				stackedSizes.add(new NOAOR(list.size(),reference,key));
				nonStackedOnReference.putIfAbsent(reference, list.size());
				for(Alignment al : list){
					String name = al.getReadName()+al.getSequence();
					if(useAllAlignments) {
						if(map.containsKey(name)){
							ArrayList<Alignment> als = map.get(name);
							als.add(al);
							map.replace(name, als);
						}else{
							ArrayList<Alignment> als = new ArrayList<Alignment>();
							als.add(al);
							map.put(al.getReadName()+al.getSequence(), als);
						}
					}else {
						if(al.isTopAlignment()) {
							ArrayList<Alignment> als = new ArrayList<Alignment>();
							als.add(al);
							map.putIfAbsent(name, als);
						}
					}
				}
			}else{
				
				GetStackedReads reads = new GetStackedReads(rMap.get(reference));
				reads.calculateStatistics();
				if(!reads.wasTurnedOn())
					turnedOn = false;
				ArrayList<Alignment> list =  reads.getNonStacked();
				stackedSizes.add(new NOAOR(list.size(),reference,key));
				for(Alignment al : list){
					if(useAllAlignments) {
						String name = al.getReadName()+al.getSequence();
						if(map.containsKey(name)){
							ArrayList<Alignment> als = map.get(name);
							als.add(al);
							map.replace(name, als);
						}else{
							ArrayList<Alignment> als = new ArrayList<Alignment>();
							list.add(al);
							map.put(name, als);
						}
					}else {
						if(al.isTopAlignment()){
							ArrayList<Alignment> als = new ArrayList<Alignment>();
							als.add(al);
							map.putIfAbsent(al.getReadName()+al.getSequence(), als);
						}
					}
				}
			}
		}
	}
	resultsMap=map;
}


// process compostion map
public void process(){
	HashMap<Integer,HashMap<String,ArrayList<Alignment>>> cMap = compositionMap;
	int max=0; 
	int maxKey=0; 
	String maxReference="";
	for(Integer id: cMap.keySet()){
		HashMap<String,ArrayList<Alignment>> map = cMap.get(id);
		for(String key : map.keySet()){
			if(map.get(key).size()>max){
				maxKey = id;
				max= map.get(key).size();
				maxReference=key;
			}
		}
	}
	//System.out.println(maxKey);
	setMaxID(maxKey);
	setMaxReference(maxReference);
	if(behave!=Filter.SRNA) {
		if(!turnOffDeDupping) {
			markAllDuplicates();
		}
	 }else {//if Filter.SNRNA we can stop destacking
		 turnOffDestacking=true;
	 }	
  }

//calculate and get all top references
public String getTopTenReferences(){
	stackedSizes.sort(new NOAORComparator());
	String line =getName(maxID)+";_"+maxReference+";_TOPREFPERCREADS100";
	if(stackedSizes.size()>=1){
		int i=0;
		for(NOAOR n : stackedSizes){
			if(n.getTaxID()!=maxID && n.getReference() != maxReference){
				int  percentage =(int) (((double) n.getSize()/(double) maxSize)*100);
				if(percentage <=9){
					line += "\t"+getName(n.getTaxID())+";_"+n.getReference()+";_TOPREFPERCREADS00"+percentage;
				}else if(percentage<100){
					line += "\t"+getName(n.getTaxID())+";_"+n.getReference()+";_TOPREFPERCREADS0"+percentage;
				}else{
					line += "\t"+getName(n.getTaxID())+";_"+n.getReference()+";_TOPREFPERCREADS100";
				}	
				i++;
			}
			if(i == 9)
				break;
		}
		if(i<9){
			while(i<9){
				line += "\tNA";
				i++;
			}
		}
	}
	return line;
}
String getName(int taxId){
	String name;
	if(mapReader.getNcbiIdToNameMap().get(taxId) != null)
		name = mapReader.getNcbiIdToNameMap().get(taxId).replace(' ', '_');
	else if(taxId == 0)
		name="NA";
	else
		name = "unassignedName";
	return name;
}
}//class
