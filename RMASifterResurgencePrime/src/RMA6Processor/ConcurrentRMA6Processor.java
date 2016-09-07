package RMA6Processor;

import java.util.List;
import java.util.concurrent.Callable;

import NCBI_MapReader.NCBI_MapReader;
import NCBI_MapReader.NCBI_TreeReader;
import behaviour.Filter;
import behaviour.Taxas;

// should be a callable version to process a RMA6 File and return 
public class ConcurrentRMA6Processor implements Callable<RMA6Processor>{
	private String outDir;
	private String fileName;
	private String inDir;
	private NCBI_MapReader mapReader;
	private NCBI_TreeReader treeReader;
	private List<Integer> taxIDs;
	private double topPercent;
	private Filter behave;
	private int maxLength;
	private Taxas t;
	private boolean readInf;
	public ConcurrentRMA6Processor(String inDir, String fileName, String outDir, NCBI_MapReader mapReader, NCBI_TreeReader treeReader,List<Integer>taxIDs,
			double topPercent, int i, Filter b, Taxas t, boolean read) {
		this.inDir = inDir;
		this.outDir = outDir;
		this.fileName = fileName;
		this.mapReader = mapReader;
		this.treeReader = treeReader;
		this.taxIDs = taxIDs;
		this.topPercent=topPercent;
		this.behave = b;
		this.maxLength = i;
		this.t = t;
		this.readInf = read;
	}
	@Override
	public RMA6Processor call(){
		RMA6Processor processor = new RMA6Processor(inDir, fileName, outDir, mapReader, treeReader,maxLength ,behave, t); // should be implemented as callable 
    	processor.process(taxIDs, topPercent, readInf);// loop through file
		return processor;
	}

}
