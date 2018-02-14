package RMA6Processor;

import java.util.ArrayList;
import java.util.concurrent.Callable;
import java.util.logging.Logger;
import NCBI_MapReader.NCBI_MapReader;
import NCBI_MapReader.NCBI_TreeReader;
import utility.InputParameterProcessor;

/**
 * Concurrent class of RMA6 Processor that allows one RMA6Processor to run per thread
 * which speeds up the process of analyzing multiple files. This seems to be the best strategy for
 * parallelization. Trying to multi-thread the process of analyzing a node severly slows down Maltextract
 * So up to date this is the fasteest way to do it
 * @author huebler
 *
 */
public class ConcurrentRMA6Processor implements Callable<RMA6Processor>{
	/**
	 * @param String indDir, String fileName, String outDir, NCBI_MapReader, NCBI_TreeReader, Set<Integer>, Filter behave,
	 * Taxas taxas, boolean verbose, int maxLength, Logger log, Logger warning
	 * @return RMA6Processor
	 * @throws none thrown all caughteName, Taxas enum, List<Intefer> UserIds, NCBI_TreeReader reader,Logger log, Logger warning
	 */
	private InputParameterProcessor inputParameterProcessor;
	private ArrayList<Integer> taxIDs;
	private NCBI_MapReader mapReader;
	private NCBI_TreeReader treeReader;
	private Logger log;
	private Logger warning;
	private String inDir;
	private String fileName;
	public ConcurrentRMA6Processor(InputParameterProcessor inputParameterProcessor,String inDir, String fileName, ArrayList<Integer> taxIDs, NCBI_MapReader mapReader,NCBI_TreeReader treeReader,
			Logger log, Logger warning) {
		
		this.inputParameterProcessor=inputParameterProcessor;
		this.taxIDs= taxIDs;
		this.mapReader=mapReader;
		this.treeReader=treeReader;
		this.log = log;
		this.warning = warning;
		this.inDir = inDir;
		this.fileName = fileName;
	}
	@Override
	public RMA6Processor call(){
		RMA6Processor processor = new RMA6Processor(inputParameterProcessor, inDir, fileName, taxIDs, mapReader, treeReader,log, warning); // should be implemented as callable 
    	processor.process();
		return processor;
	}

}
