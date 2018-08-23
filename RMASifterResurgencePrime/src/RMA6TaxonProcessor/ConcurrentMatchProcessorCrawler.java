package RMA6TaxonProcessor;

import java.util.ArrayList;
import java.util.concurrent.Callable;
import java.util.logging.Logger;

import NCBI_MapReader.NCBI_MapReader;
import RMAAlignment.Alignment;

public class ConcurrentMatchProcessorCrawler implements Callable<MatchProcessorCrawler> {
	/**
	 * Is used to somewhat parallelize RMA6Crawler. When RMA6 Reads are sorted to the correct 
	 * species one node per thread can be processed with this class. some values are initialized 
	 * at default value. 
	 * @author huebler
	 * @params ConcurrentLinkedDeque<Alignment> concurrentLinkedDeque,int id,Logger log, Logger warning,
	 * boolean wantReads, NCBI_MapReader mapReader
	 */
	private int taxID;
	private double top = 0.01;

	private Logger log;
	private Logger warning;
	private ArrayList<Alignment> clD; 
	private boolean wantReads;
	private boolean singleStranded;
	
	private NCBI_MapReader mapReader;
	public ConcurrentMatchProcessorCrawler(ArrayList<Alignment> concurrentLinkedDeque,int id,Logger log, Logger warning,boolean wantReads,
			NCBI_MapReader mapReader, boolean singleStranded){
		this.log = log;
		this.warning= warning;
		this.wantReads = wantReads;
		this.taxID = id;
		this.mapReader = mapReader;
		clD=concurrentLinkedDeque;
		this.singleStranded = singleStranded;
	}
	@Override
	public MatchProcessorCrawler call() throws Exception {
		System.out.println("size list " + clD.size());
		MatchProcessorCrawler matchProcessorCrawler = new MatchProcessorCrawler(taxID,top,mapReader,false,log,warning,wantReads,0.01,0,false,true,true,behaviour.Filter.CRAWL,false,singleStranded);
		matchProcessorCrawler.processDLQlist(clD);
		matchProcessorCrawler.process();
		return matchProcessorCrawler;
	}

}
