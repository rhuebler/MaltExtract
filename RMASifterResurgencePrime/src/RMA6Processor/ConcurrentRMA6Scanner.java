package RMA6Processor;

import java.util.List;
import java.util.concurrent.Callable;

import NCBI_MapReader.NCBI_TreeReader;
import behaviour.Taxas;

public class ConcurrentRMA6Scanner implements Callable<RMA6Scanner>{
	private String inDir;
	private String fileName;
	private List<Integer> TaxIDs;
	private Taxas tax;
	private NCBI_TreeReader tReader;
	public ConcurrentRMA6Scanner(String inDir, String fileName,Taxas tax,List<Integer> Ids,NCBI_TreeReader tReader){
		this.inDir = inDir;
		this.fileName = fileName;
		this.TaxIDs =Ids;
		this.tax = tax;
		this.tReader = tReader;
	}
	@Override
	public RMA6Scanner call(){
		RMA6Scanner scanner = new RMA6Scanner(inDir, fileName,tax,TaxIDs, tReader); // should be implemented as callable 
		return scanner;
	}
}
