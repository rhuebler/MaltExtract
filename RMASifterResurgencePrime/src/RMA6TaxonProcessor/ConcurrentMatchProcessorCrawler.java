package RMA6TaxonProcessor;

import java.util.concurrent.Callable;

public class ConcurrentMatchProcessorCrawler implements Callable<MatchProcessorCrawler> {
	MatchProcessorCrawler matchProcessorCrawler; 
	public ConcurrentMatchProcessorCrawler(MatchProcessorCrawler mpc){
		matchProcessorCrawler=mpc;
	}
	@Override
	public MatchProcessorCrawler call() throws Exception {
		matchProcessorCrawler.process();
		return matchProcessorCrawler;
	}

}
