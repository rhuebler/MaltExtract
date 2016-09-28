package RMA6TaxonProcessor;

import java.util.concurrent.Callable;

public class ConcurrentRMA6TaxonProcessor implements Callable<RMA6TaxonProcessor> {
	private RMA6TaxonProcessor taxProcessor;
	private String inDir;
	private String fileName;
	private double topPercent;
	private int maxLength;
	public ConcurrentRMA6TaxonProcessor(RMA6TaxonProcessor taxProcessor, String s1, String s2, double tp, int mL){
	this.taxProcessor = taxProcessor;
	this.inDir = s1;
	this.fileName = s2;
	this.topPercent = tp;
	this.maxLength = mL;
	}
	
	@Override
	public RMA6TaxonProcessor call(){
		try{
		taxProcessor.process(inDir, fileName, topPercent, maxLength);
		}catch(Exception e){
			e.printStackTrace();
		}
		return this.taxProcessor;
	}
}