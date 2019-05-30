package NCBI_MapReader;

import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.net.URL;
import java.net.URLConnection;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.zip.ZipEntry;
import java.util.zip.ZipInputStream;

public class testy {
	static HashMap<Integer,ArrayList<String>> IDsToNames = new HashMap<Integer,ArrayList<String>>();
	public static void processNCBIZipFile() {
		
		String location = "ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdmp.zip";
		try{
			byte[] array = null;
			URLConnection conn = new URL(location).openConnection();
			 conn.setConnectTimeout(90*1000);
			 conn.setReadTimeout(90*1000);
			   try (InputStream in = conn.getInputStream()) {
				   ZipInputStream zipStream = new ZipInputStream(in);
				   ZipEntry zipEntry = zipStream.getNextEntry();
			       while (zipEntry != null) {
			    	   if(zipEntry.getName().equals("names.dmp")) {
			    		   byte[] btoRead = new byte[1024];
				    	   ByteArrayOutputStream bout = new ByteArrayOutputStream();; //<- I don't want this!
				            int len = 0;
				            while ((len = zipStream.read(btoRead)) != -1) {
				                bout.write(btoRead, 0, len);
				            }
				            bout.close();
				            array =  bout.toByteArray();
			    	   }
			    	
			            zipEntry = zipStream.getNextEntry();
			       }
				   zipStream.close();
				  
				 String lines = new String(array); 
				 
				 for(String line:lines.split("\\n")) {
					String[] parts = line.split("\\t");
					//System.out.println(parts[0]+"\t"+parts[2]);
					if(!IDsToNames.containsKey(Integer.parseInt(parts[0]))) {
						ArrayList<String> entries = new ArrayList<String>();
						entries.add(parts[2]);
						IDsToNames.put(Integer.parseInt(parts[0]), entries);
					}	
					else {
						ArrayList<String> entries = IDsToNames.get(Integer.parseInt(parts[0]));
						entries.add(parts[2]);
						IDsToNames.replace(Integer.parseInt(parts[0]),entries);
					}
				 }
				 System.out.println("Done loading names");
			   }catch(Exception e) {
				   e.printStackTrace();
			    }
		}catch(IOException io) {
			io.printStackTrace();
		}
		
	}
 public static void main(String[] args) {
		NCBI_MapReader mapReader;
		NCBI_TreeReader treeReader;
		mapReader = new NCBI_MapReader();
		treeReader = new NCBI_TreeReader();
		processNCBIZipFile();
		System.out.println("Done loading tree and ncbi map");
		for(int ID: treeReader.getLeavesIds()) {
			//System.out.println(ID);
			if(IDsToNames.containsKey(ID)) {
				ArrayList<String> entries = IDsToNames.get(ID);
				if(entries.size()>1) {
					System.out.print(ID +"	"+mapReader.getNcbiIdToNameMap().get(ID));
					for(String entry : entries) {
						System.out.print("	"+entry);
					}
					System.out.println("");
				}else {
					System.out.println(ID +"	"+mapReader.getNcbiIdToNameMap().get(ID)+"	"+entries.get(0));
				}
			}else {
				System.out.println(ID+"	"+mapReader.getNcbiIdToNameMap().get(ID));
			}
		}
 	}
}
