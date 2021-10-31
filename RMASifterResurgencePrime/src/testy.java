import java.util.ArrayList;
import java.util.Set;
import java.util.logging.Logger;

import NCBI_MapReader.NCBI_MapReader;
import NCBI_MapReader.NCBI_TreeReader;
import RMAExtractor.RMAExtractor;
import megan.rma6.ClassificationBlockRMA6;
import megan.rma6.RMA6File;

public class testy {
	private static Set<Integer> keySet;
	private static Set<Integer> allKeys;
	private static Integer readCount;
	public static void main(String[] args) {
		try{

		Logger logger = Logger.getLogger(testy.class.getName());
		Logger error= Logger.getLogger("Error");
		NCBI_TreeReader treeReader = new NCBI_TreeReader(logger, error);
		NCBI_MapReader	mapReader = new NCBI_MapReader(logger, error);	
		String fileName = "/Users/huebler/Desktop/Edaphobacter_sp._12200R-103_strain_12200R-103_ASM1009302v1_2703788.fna.gz.simulated.rma6";
		String parts[] = fileName.split("_");
    	int taxID= Integer.parseInt(parts[(parts.length-1)].split("\\.")[0]);
    	System.out.println(taxID);
		RMA6File rma6File = new RMA6File(fileName, "r");
		Long location = rma6File.getFooterSectionRMA6().getStartClassification("Taxonomy");

	    if (location != null) {
	        ClassificationBlockRMA6 cl = new ClassificationBlockRMA6("Taxonomy");
	        cl.read(location, rma6File.getReader());
	      
	        keySet = cl.getKeySet();
	        //this.allKeys = cl.getKeySet();
	        readCount = (int) rma6File.getFooterSectionRMA6().getNumberOfReads();// read in all counts
			if(mapReader.getNcbiIdToNameMap().containsKey(taxID)) {
	        ArrayList<Integer> onPathIDs = treeReader.getTaxonomicPath(taxID, keySet);
	        for(int id:onPathIDs) {
	        	System.out.println(id);

	        }
			}
	}
	}catch(Exception e) {e.printStackTrace();}
	}
}
