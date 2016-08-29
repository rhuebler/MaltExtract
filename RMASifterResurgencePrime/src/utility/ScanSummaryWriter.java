package utility;

import java.io.IOException;
import java.nio.charset.Charset;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.List;
import java.util.Set;

import NCBI_MapReader.NCBI_MapReader;
import RMA6Processor.RMA6Scanner;
/**
 * Write a summary File to quickly summarize which IDs have assigned nodes within a RMA6 file and how many
 * seems to be very fast even on one core longest part is loading the NCBI.tre an NCBI.map for IDs and taxonomy 
 * @author huebler
 *
 */
public class ScanSummaryWriter {
	private List<RMA6Scanner> list;
	private Set<Integer> keys;
	private NCBI_MapReader reader;
	public ScanSummaryWriter(List<RMA6Scanner> scannerList,Set<Integer> allKeys, NCBI_MapReader reader){
		this.list = scannerList;
		this.keys = allKeys;
		this.reader = reader;
	}
	public void write(String outDir){
		boolean first = true;
		List<String> summary = new ArrayList<String>();
		String header = "Node_Name";
		for(RMA6Scanner scan : list){
			header += "\t" + scan.getFileName();
			if(first){
				for(int key : keys){
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
				for(int key : keys){
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
		summary.add(0,header);
		try{
		System.out.println("Writing Scan_Summary txt File");
		Path file = Paths.get(outDir+"ScanSummary"+".txt");
		Files.write(file, summary, Charset.forName("UTF-8"));
		System.out.println("Scan_Summary Done!");
		}catch(IOException io){
			io.printStackTrace();
		}
	}
	
	
}
