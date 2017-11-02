package utility;

import java.util.concurrent.Callable;
import java.util.logging.Logger;
import NCBI_MapReader.NCBI_MapReader;


public class ConcurrentNodeMatchSorter implements Callable<NodeMatchSorter> {
	private String filePath="";
	int id =0;
	Logger log;
	Logger warning;
	boolean wantReads;
	double topPercent;
	String speciesName;
	NCBI_MapReader mapReader;
	public ConcurrentNodeMatchSorter(String filePath, int id,Logger log, Logger warning,boolean wantReads,
	String speciesName,NCBI_MapReader mapReader){
		this.filePath = filePath;
		this.id = id;
		this.warning = warning;
		this.speciesName = speciesName;
		this.mapReader=mapReader;
		this.log = log;
		this.wantReads =wantReads;
	}
	
	public NodeMatchSorter call() {
		NodeMatchSorter nodeMatchSorter = new NodeMatchSorter(filePath, id, log, warning, wantReads, speciesName, mapReader);
		try{
			nodeMatchSorter.processNode();
			
		}catch(Exception e){
			e.printStackTrace();
		}
		return nodeMatchSorter;
	}
}
