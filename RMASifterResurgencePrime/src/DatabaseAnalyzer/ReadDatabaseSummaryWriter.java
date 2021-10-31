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
				header+="\tOnPathPercentageStrict\tOffPathPercentageStrict\tOnPathPercentageRelaxed\tOffPathPercentagRelaxed\tisMonoCladic\tTotalCount\t"
						+ "Node_01\tNode_02\tNode_03\tNode_04\tNode_05\tNode_06\tNode_07\tNode_08\tNode_09\tNode_10";
				for(String name : fileNames) {
					if(map.containsKey(name)) {
						try {
							ReadDatabaseAnalyzer analyzer=map.get(name).get();
							output.add(new File(name).getName()+"\t"+analyzer.getOnPathStrict()+"\t"+analyzer.getOffPathStrict()+"\t"+analyzer.getOnPathRelaxed()+"\t"
							+analyzer.getOffPathRelaxed()+"\t"+analyzer.isMonoCladic()+"\t"+analyzer.getTotalCount()+"\t"+analyzer.getOffPathNodes());
						} catch (InterruptedException | NullPointerException | ExecutionException e) {
							// TODO Auto-generated catch block
							output.add(new File(name).getName()+"\tNA\tNA\tNA\tNA\tNA\tNA");
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
