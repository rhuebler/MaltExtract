package DatabaseAnalyzer;

import java.io.File;
import java.io.IOException;
import java.nio.charset.Charset;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.Future;
import java.util.logging.Level;
import java.util.logging.Logger;


public class ReadDatabaseSummaryWriter {

	/**
	 * @param List<Future<RMA6Scanner>>, NCBI_MapReader, Set<Integer>, Logger,Warning
	 * @throws none thrown all caught
	 */
	private HashMap<String,Future<ReadDatabaseAnalyzer>> map;

	private List<String> summary;
	private Logger warning;
	private  ArrayList<String> fileNames;
	private DatabaseAnalysisMode mode;
	public ReadDatabaseSummaryWriter(HashMap<String,Future<ReadDatabaseAnalyzer>> databaseMap,Logger warning, ArrayList<String> fileNames, DatabaseAnalysisMode mode){
		this.map = databaseMap;
		this.warning = warning;
		this.fileNames = fileNames;
		this.mode = mode;
		prepareOutput();
	}
	private void prepareOutput(){
		String header = "FileName";
		ArrayList<String> output= new ArrayList<String>();
		fileNames.sort(null);
		switch(mode) {
			case LIST:{
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
			break;
			}
			case ONPATH:{
				header+="\tOnPathPercentage\tOffPathPercentage";
				for(String name : fileNames) {
					if(map.containsKey(name)) {
						try {
							ReadDatabaseAnalyzer analyzer=map.get(name).get();
							
							output.add(new File(name).getName()+"\t"+analyzer.getOnPath()+"\t"+analyzer.getOffPath());
						} catch (InterruptedException | ExecutionException e) {
							// TODO Auto-generated catch block
							e.printStackTrace();
						}
					}		
				}
				break;
			}
		}
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
