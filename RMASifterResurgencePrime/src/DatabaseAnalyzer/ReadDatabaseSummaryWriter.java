package DatabaseAnalyzer;

import java.io.IOException;
import java.nio.charset.Charset;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Set;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.Future;
import java.util.logging.Level;
import java.util.logging.Logger;

import NCBI_MapReader.NCBI_MapReader;

public class ReadDatabaseSummaryWriter {

	/**
	 * @param List<Future<RMA6Scanner>>, NCBI_MapReader, Set<Integer>, Logger,Warning
	 * @throws none thrown all caught
	 */
	private HashMap<String,Future<ReadDatabaseAnalyzer>> map;

	private List<String> summary;
	private Logger warning;
	private  ArrayList<String> fileNames;
	public ReadDatabaseSummaryWriter(HashMap<String,Future<ReadDatabaseAnalyzer>> databaseMap,Logger warning, ArrayList<String> fileNames){
		this.map = databaseMap;
		this.warning = warning;
		this.fileNames = fileNames;
		prepareOutput();
	}
	private void prepareOutput(){
		String header = "FileName\t";
		ArrayList<String> output= new ArrayList<String>();
		fileNames.sort(null);
		for(int i =1; i<=10; i++) {
			header+="\tNode_"+i;
			}
		for(String name : fileNames) {
			if(map.containsKey(name)) {
				try {
					output.add(map.get(name).get().getOutput());
				} catch (InterruptedException | ExecutionException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
			}else {
				int i=0;
				String line = "name\t";
				while(i<10) {
					line+="NA;NA\t";
				}
				output.add(line);
			}
		}
		output.sort(null);
		output.add(0,header);
		summary = output;
	}
		public void write(String outDir){
			try{
				Path file = Paths.get(outDir+"AssingedNodes"+".txt");
				Files.write(file, summary, Charset.forName("UTF-8"));
			}catch(IOException io){
				warning.log(Level.SEVERE, "Error", io);
			}
		}
	}
