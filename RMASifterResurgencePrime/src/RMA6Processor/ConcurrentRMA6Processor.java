package RMA6Processor;

import java.util.List;
import java.util.concurrent.Callable;

import NCBI_MapReader.NCBI_MapReader;
import NCBI_MapReader.NCBI_TreeReader;
import behaviour.Behaviour;

// should be a callable version to process a RMA6 File and return 
public class ConcurrentRMA6Processor implements Callable<RMA6Processor>{
	private String outDir;
	private String fileName;
	private String inDir;
	private NCBI_MapReader mapReader;
	private NCBI_TreeReader treeReader;
	private List<Integer> taxIDs;
	private double topPercent;
	private Behaviour behave;
	public ConcurrentRMA6Processor(String inDir, String fileName, String outDir, NCBI_MapReader mapReader, NCBI_TreeReader treeReader,List<Integer>taxIDs, double topPercent, Behaviour b) {
		this.inDir = inDir;
		this.outDir = outDir;
		this.fileName = fileName;
		this.mapReader = mapReader;
		this.treeReader = treeReader;
		this.taxIDs = taxIDs;
		this.topPercent=topPercent;
		this.behave = b;
	}
	@Override
	public RMA6Processor call() throws Exception {
		RMA6Processor processor = new RMA6Processor(inDir, fileName, outDir, mapReader, treeReader, behave); // should be implemented as callable 
    	processor.process(taxIDs, topPercent);// loop through file
		return processor;
	}

}
