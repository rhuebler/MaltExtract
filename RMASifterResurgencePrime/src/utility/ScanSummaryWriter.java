package utility;

import java.io.IOException;
import java.nio.charset.Charset;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.Future;
import java.util.logging.Level;
import java.util.logging.Logger;

import NCBI_MapReader.NCBI_MapReader;
import RMA6Processor.RMA6Scanner;
/**
 * Write a summary File to quickly summarize which IDs have assigned nodes within a RMA6 file and how many
 * seems to be very fast even on one core longest part is loading the NCBI.tre an NCBI.map for IDs and taxonomy 
 * @author huebler
 *
 */
public class ScanSummaryWriter {
	/**
	 * @param List<Future<RMA6Scanner>>, NCBI_MapReader, Set<Integer>, Logger,Warning
	 * @throws none thrown all caught
	 */
	private List<Future<RMA6Scanner>> list;
	private NCBI_MapReader reader;
	private Set<Integer> keySet;
	private List<String> summary;
	private Logger warning;
	public ScanSummaryWriter(List<Future<RMA6Scanner>> scannerList, NCBI_MapReader reader,Logger warning){
		this.list = scannerList;
		this.reader = reader;
		this.warning = warning;
		setProcessedIds();
		prepareOutput();
	}
	void setProcessedIds(){
		Set<Integer> keys = new HashSet<>();
		try{
		for(Future<RMA6Scanner> future : list)
			keys.addAll(future.get().getKeySet());
		this.keySet = keys;
		}catch(InterruptedException ie){
			warning.log(Level.SEVERE, "Error", ie);
		}catch(ExecutionException ee){
			warning.log(Level.SEVERE, "Error", ee);
		}
	}
	private void prepareOutput(){
		List<String> summary = new ArrayList<String>();
		try{
			boolean first = true;
			String header = "Node_Name";
			for(Future<RMA6Scanner> future : list){
				RMA6Scanner scan = future.get();
				header += "\t" + scan.getFileName();
				if(first){
					for(int key : keySet){
						String s;
						if( reader.getNcbiIdToNameMap().get(key) != null){
							s = reader.getNcbiIdToNameMap().get(key).replace(' ', '_');
						}else{
							s = "unasigned";
						}
						if(scan.getKeySet().contains(key)){
							s += "\t"+scan.getAssignmentMap().get(key);
						}else{
							s += "\t0";
						}
						summary.add(s);
					}
					first=false;
				}else{
					int i = 0;
					for(int key : keySet){
						String s = summary.get(i);
						if(scan.getKeySet().contains(key)){
							s += "\t"+scan.getAssignmentMap().get(key);
						}else{
							s += "\t0";
						}
						summary.set(i, s);
						i++;
					}
				}
			}
			summary.sort(null);
			summary.add(0,header);
			this.summary =  summary;
		}catch(InterruptedException ie){
			warning.log(Level.SEVERE, "Error", ie);
		}catch(ExecutionException ee){
			warning.log(Level.SEVERE, "Error", ee);
		}
	}
		public void write(String outDir){
			try{
				Path file = Paths.get(outDir+"ScanSummary"+".txt");
				Files.write(file, summary, Charset.forName("UTF-8"));
			}catch(IOException io){
				warning.log(Level.SEVERE, "Error", io);
			}
		}
	}
	

