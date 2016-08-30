package RMA6Processor;

import java.util.concurrent.Callable;

public class ConcurrentRMA6Scanner implements Callable<RMA6Scanner>{
	private String inDir;
	private String fileName;
	public ConcurrentRMA6Scanner(String inDir, String fileName){
		this.inDir = inDir;
		this.fileName = fileName;
	}
	@Override
	public RMA6Scanner call(){
		RMA6Scanner scanner = new RMA6Scanner(inDir, fileName); // should be implemented as callable 
		return scanner;
	}
}
