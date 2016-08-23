package RMA6Processor;

import java.io.IOException;

import megan.rma6.RMA6Connector;

/**
 * Class returns the keyset of a RMA6File
 * which means it returns the names of all Nodes that have an assigned Name
 * @author huebler
 *
 */
public class RMA6Scanner {
	private String inDir;
	private String fileName;
	public RMA6Scanner(String inDir, String name){
	this.inDir = inDir;
	this.fileName =  name;
	process();
	}
	private void process(){
		try{
			RMA6Connector fileCon = new RMA6Connector(inDir+fileName);
		}catch(IOException io){
			io.printStackTrace();
		}
	}
}
