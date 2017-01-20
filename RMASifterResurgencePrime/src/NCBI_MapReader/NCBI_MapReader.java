package NCBI_MapReader;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.HashMap;
import java.util.Scanner;
/**
 * 
 * @author huebler
 * process the NCBI.map file 
 * into Map of ID to name
 * and a map of name to ID
 */

import utility.ResourceFinder;
public class NCBI_MapReader {
	/**
	 * @param String directory/to/file
	 * @throws none thrown all caught
	 * @return return two hashmaps with IDs to name and Name to ID
	 */
	private String mapName = "";
	 HashMap<String,Integer> ncbiNameToId;
	 HashMap<Integer,String> ncbiIdToName;
	public NCBI_MapReader(){
		ResourceFinder resources = new ResourceFinder();
		mapName = resources.getPath("ncbi.map");
		processNcbiMap(mapName);
	}
	public NCBI_MapReader(String path){
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
	}