package DatabaseAnalyzer;
import java.util.concurrent.Callable;
import java.util.logging.Logger;
import NCBI_MapReader.NCBI_MapReader;

public class ConcurrentReadDatabaseAnalyzer implements Callable<ReadDatabaseAnalyzer> {
	private String inDir;
	private String fileName ;
	private Logger log ;
	private Logger warning ;
	private NCBI_MapReader reader;
	public ConcurrentReadDatabaseAnalyzer(String inDir, String name,
			Logger log, Logger warning, NCBI_MapReader reader){
		this.inDir = inDir;
		this.fileName =  name;
		this.log = log;
		this.warning = warning;
		this.reader = reader;
		
	}
	@Override
	public ReadDatabaseAnalyzer call() throws Exception {
		ReadDatabaseAnalyzer analyzer = new ReadDatabaseAnalyzer(inDir, fileName, log, warning, reader);
		
		// TODO Auto-generated method stub
		return analyzer;
	}

}
