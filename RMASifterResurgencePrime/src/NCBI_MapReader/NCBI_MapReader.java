package NCBI_MapReader;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.net.URL;
import java.net.URLConnection;
import java.util.HashMap;
import java.util.Scanner;
/**
 * 
 * @author huebler
 * process the NCBI.map file 
 * into Map of ID to name
 * and a map of name to ID
 */

public class NCBI_MapReader {
	/**
	 * Is used to read the NCBI map fie and generate two hashmaps that 
	 * map the taxonomy ID to the taxonomy name and the name to the taxonomy ID
	 * which is necessary to interface human understandable species names to the 
	 * taxonomy IDs used internally to exit RMA6 files
	 * @param String directory/to/file
	 * @throws none thrown all caught
	 * @return return two hashmaps with IDs to name and Name to ID
	 */
	private String mapName = "";//TODO download file
	 HashMap<String,Integer> ncbiNameToId;
	 HashMap<Integer,String> ncbiIdToName;
	public NCBI_MapReader(){// try to locate resources
		processFromWeb();
	}
	public NCBI_MapReader(String path){// use provided path
		if(!path.endsWith("/"))
			path+="/";
		this.mapName = path + "ncbi.map";
		processNcbiMap(mapName);
	}
	public HashMap<String,Integer>  getNcbiNameToIdMap(){
		return this.ncbiNameToId;
	}
	
	public HashMap<Integer,String>  getNcbiIdToNameMap(){
		return this.ncbiIdToName;
	}
	
	public void  setNcbiNameToIdMap(HashMap<String,Integer> map){
		this.ncbiNameToId = map;
	}
	
	private void  setNcbiIdToNameMap(HashMap<Integer,String> map){
		this.ncbiIdToName = map;
	}
	
	public void processNcbiMap(String fileName) {
		try(Scanner in = new Scanner(new File(fileName))){
			HashMap<String,Integer> ncbiNameMap = new HashMap<String,Integer>();
			HashMap<Integer, String> ncbiIDMap = new HashMap<Integer, String>();
			while (in.hasNext()) { // iterates each line in the file
				String line = in.nextLine();
				String[] frags = line.split("\\t");
				ncbiNameMap.put(frags[1], Integer.parseInt(frags[0]));
				ncbiIDMap.put(Integer.parseInt(frags[0]), frags[1]);
				// do something with line
				}
			in.close();
			setNcbiNameToIdMap(ncbiNameMap);
			setNcbiIdToNameMap(ncbiIDMap);
			}catch(FileNotFoundException ie){
				ie.printStackTrace();
			}
		}
		public void processFromWeb() {
			try{
				HashMap<String,Integer> ncbiNameMap = new HashMap<String,Integer>();
				HashMap<Integer, String> ncbiIDMap = new HashMap<Integer, String>();
				//https://github.com/danielhuson/megan-ce/blob/master/resources/files/ncbi.map
				String location = "https://raw.githubusercontent.com/danielhuson/megan-ce/master/resources/files/ncbi.map";
				URLConnection conn = new URL(location).openConnection();
				 conn.setConnectTimeout(90*1000);
				 conn.setReadTimeout(90*1000);
				   try (InputStream in = conn.getInputStream()) {
					   InputStreamReader reader = new  InputStreamReader(in);
					   BufferedReader buffered = new BufferedReader(reader);
					   String line;
					   while((line = buffered.readLine())!=null) {
						   System.out.println(line);
							//String[] frags = line.toString().split("\\t");
//							ncbiNameMap.put(frags[1], Integer.parseInt(frags[0]));
//							ncbiIDMap.put(Integer.parseInt(frags[0]), frags[1]);
							// do something with line
				       }
					 buffered.close();
					setNcbiNameToIdMap(ncbiNameMap);
					setNcbiIdToNameMap(ncbiIDMap);
				   }catch(Exception e) {
					   e.printStackTrace();
				    }
			}catch(IOException io) {
				io.printStackTrace();
			}
		}
	}