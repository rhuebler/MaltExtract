package RMA6TaxonProcessor;

import java.util.concurrent.Callable;

public class ConcurrentNodeProcessor implements Callable<NodeProcessor> {
	private NodeProcessor nodeProcessor;
	private String inDir;
	private String fileName;
	private double topPercent;
	private int maxLength;
	public ConcurrentNodeProcessor(NodeProcessor nodeProcessor, String s1, String s2, double tp, int mL){
	this.nodeProcessor = nodeProcessor;
	this.inDir = s1;
	this.fileName = s2;
	this.topPercent = tp;
	this.maxLength = mL;
	}
	
	@Override
	public NodeProcessor call(){
		try{
		nodeProcessor.process(inDir, fileName, topPercent, maxLength);
		}catch(Exception e){
			e.printStackTrace();
		}
		return this.nodeProcessor;
	}
}