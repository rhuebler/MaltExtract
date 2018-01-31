package RMA6TaxonProcessor;

import java.util.concurrent.Callable;
import java.util.concurrent.ConcurrentLinkedDeque;
import java.util.logging.Logger;

import NCBI_MapReader.NCBI_MapReader;
import RMAAlignment.Alignment;

public class ConcurrentMatchProcessorCrawler implements Callable<MatchProcessorCrawler> {
	private int taxID;
	private double top = 0.01;

	private Logger log;
	private Logger warning;
	private ConcurrentLinkedDeque<Alignment> clD; 
	private boolean wantReads;
	
	private NCBI_MapReader mapReader;
	public ConcurrentMatchProcessorCrawler(ConcurrentLinkedDeque<Alignment> concurrentLinkedDeque,int id,Logger log, Logger warning,boolean wantReads,
			NCBI_MapReader mapReader){
		this.log = log;
		this.warning= warning;
		this.wantReads = wantReads;
		this.taxID = id;
		this.mapReader = mapReader;
		clD=concurrentLinkedDeque;
	}
	@Override
	public MatchProcessorCrawler call() throws Exception {
		MatchProcessorCrawler matchProcessorCrawler = new MatchProcessorCrawler(taxID,top,mapReader,false,log,warning,wantReads,0.01,0,false,true,false,behaviour.Filter.CRAWL);
		matchProcessorCrawler.processDLQlist(clD);
		matchProcessorCrawler.process();
		return matchProcessorCrawler;
	}

}
