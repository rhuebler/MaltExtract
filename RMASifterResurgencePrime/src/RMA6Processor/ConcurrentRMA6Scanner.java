package RMA6Processor;

import java.util.List;
import java.util.concurrent.Callable;
import java.util.logging.Logger;

import NCBI_MapReader.NCBI_TreeReader;
import behaviour.Taxas;
/**
 * Concurrent version of RMA6Scanner that allows for one to run on each thread
 * @author huebler
 *
 */
public class ConcurrentRMA6Scanner implements Callable<RMA6Scanner>{
	/**
	 * @param String inDir, String FileName, String outDir, List<Integer> taxIDs
	 * Taxas tax,TreeReader reader, Logger log, Logger warning
	 * @return ap<Integer,Integer> assignmentMap
	 * @throws none thrown all caughteName, Taxas enum, List<Intefer> UserIds, NCBI_TreeReader reader,Logger log, Logger warning
	 */
	private String inDir;
	private String fileName;
	private List<Integer> TaxIDs;
	private Taxas tax;
	private NCBI_TreeReader tReader;
	private Logger log;
	private Logger warning;
	private String outDir;
	private boolean wantMeganSummaries;
	public ConcurrentRMA6Scanner(String inDir, String fileName,Taxas tax,List<Integer> Ids,NCBI_TreeReader tReader,
			Logger log, Logger warning, String outDir, boolean wantMeganSummaries){
		this.inDir = inDir;
		this.fileName = fileName;
		this.TaxIDs =Ids;
		this.tax = tax;
		this.tReader = tReader;
		this.log = log;
		this.warning = warning;
		this.outDir = outDir;
		this.wantMeganSummaries = wantMeganSummaries;
	}
	@Override
	public RMA6Scanner call(){
		RMA6Scanner scanner = new RMA6Scanner(inDir, fileName,tax,TaxIDs, tReader, log, warning, outDir,wantMeganSummaries); // should be implemented as callable 
		return scanner;
	}
}
