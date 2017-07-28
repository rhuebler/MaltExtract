package RMA6Processor;

import java.util.List;
import java.util.concurrent.Callable;
import java.util.logging.Logger;

import NCBI_MapReader.NCBI_MapReader;
import NCBI_MapReader.NCBI_TreeReader;
import behaviour.Filter;
import behaviour.Taxas;

/**
 * Concurrent class of RMA6 Processor that allows for one RMA6Processor per thread
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
	//intilaize attributes
	private String outDir;
	private String fileName;
	private String inDir;
	private NCBI_MapReader mapReader;
	private NCBI_TreeReader treeReader;
	private List<Integer> taxIDs;
	private double topPercent;
	private Filter behave;
	private int maxLength;
	private double minPIdent;
	private Taxas t;
	private boolean verbose;
	private boolean wantReads = false;
	private Logger log;
	private Logger warning;
	private double minComplexity;
	private boolean wantAlignments;
	private  boolean wantMeganSummaries = false;
	private boolean turnOffDestacking = false;
	private boolean turnOffDeDuping = false;
	public ConcurrentRMA6Processor(String inDir, String fileName, String outDir, NCBI_MapReader mapReader, 
			NCBI_TreeReader treeReader,List<Integer>taxIDs,double topPercent, int i,double minPI, Filter b,
			Taxas t, boolean verbose, Logger log, Logger warning, boolean reads,  double minCompl, boolean alignments, boolean wantMeganSummaries, boolean turnOffDestacking, boolean dedupOff) {
		// set atrributes to value
		this.inDir = inDir;
		this.outDir = outDir;
		this.fileName = fileName;
		this.mapReader = mapReader;
		this.treeReader = treeReader;
		this.taxIDs = taxIDs;
		this.topPercent=topPercent;
		this.behave = b;
		this.maxLength = i;
		this.minPIdent = minPI;
		this.t = t;
		this.verbose = verbose;
		this.log = log;
		this.warning = warning;
		this.wantReads = reads;
		this.wantAlignments = alignments;
		this.minComplexity = minCompl;
		this.wantMeganSummaries = wantMeganSummaries;
		this.turnOffDestacking = turnOffDestacking;
		this.turnOffDeDuping = dedupOff;
	}
	@Override
	public RMA6Processor call(){
		//initialize RMA6Processor and process RMA6 file
		RMA6Processor processor = new RMA6Processor(inDir, fileName, outDir, mapReader,
				treeReader,maxLength,minPIdent ,behave, t, verbose, log, warning, wantReads, minComplexity,wantAlignments,wantMeganSummaries, turnOffDestacking,turnOffDeDuping); // should be implemented as callable 
    	processor.process(taxIDs, topPercent);// loop through file
		return processor;
	}

}
