package utility;

import java.util.concurrent.ConcurrentHashMap;
import java.util.logging.Logger;

import NCBI_MapReader.NCBI_MapReader;
import RMA6TaxonProcessor.MatchProcessorCrawler;

public class ConcurrentNodeMatchSorter implements Runnable {
	private ConcurrentHashMap<Integer, MatchProcessorCrawler> concurrentMap;
	private String filePath="";
	int id =0;
	Logger log;
	Logger warning;
	boolean wantReads;
	double topPercent;
	String speciesName;
	NCBI_MapReader mapReader;
	public ConcurrentNodeMatchSorter(ConcurrentHashMap<Integer, MatchProcessorCrawler> concurrentMap,String filePath, int id,Logger log, Logger warning,boolean wantReads,
	String speciesName,NCBI_MapReader mapReader){
		this.concurrentMap = concurrentMap;
		this.filePath = filePath;
		this.id = id;
		this.warning = warning;
		this.speciesName = speciesName;
		this.mapReader=mapReader;
		this.log = log;
		this.wantReads =wantReads;
	}
	
	@Override
	public void run() {
		try{
			NodeMatchSorter nodeMatchSorter = new NodeMatchSorter(concurrentMap, filePath, id, log, warning, wantReads, speciesName, mapReader);
			nodeMatchSorter.processNode();
		}catch(Exception e){
			e.printStackTrace();
		}
		
	}
}
