package NCBI_MapReader;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashSet;
import java.util.Scanner;
import java.util.Set;

import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;
/**
 * Process NCBI tre File with archeoptrix.jar
 * Process tree find child nodes to input taxID
 * only keep child IDs that are assigned in current RMA6 file 
 * Unfortunately does Malt only store taxIds of nodes that have at least one assigned 
 * read therefore the number of nodes cannot be reduced at every level for the risk of losing 
 * a branch. Instead a more crawling like approach is currently needed that allows to skip one taxonomic evel or two 
 * @author huebler
 *
 */

public class NCBI_TreeReader {
	/**
	 * @param String directory/to/file
	 * @throws none thrown all caught
	 * @return return phylogenetic tree object
	 */
	private String treName= "/projects1/clusterhomes/huebler/RMASifter/RMA_Extractor_Resources/ncbi.tre";//TODO how to provide System resources make relativistic
	int target;
	private HashSet<Integer> positionsToKeep = new HashSet<Integer>();
	private Phylogeny ph;
	public NCBI_TreeReader(){
		try(Scanner in = new Scanner(new File(treName))){
		String line = in.nextLine();
		in.close();
	    this.ph = Phylogeny.createInstanceFromNhxString(line);
	    }catch(IOException io){
	    	io.printStackTrace();
	    }
	}
	public NCBI_TreeReader(String path){
		this.treName = path + "ncbi.tre";
		try{
			Scanner in = new Scanner(new File(treName));
			String line = in.nextLine();
			in.close();
		    this.ph = Phylogeny.createInstanceFromNhxString(line);
		    }catch(IOException io){
		    	io.printStackTrace();
		    }
	}
	public NCBI_TreeReader(NCBI_TreeReader copyInstance){
		this.ph = copyInstance.getPhylogeny();
	}
	private Phylogeny getPhylogeny(){
		return ph;
	}
	private ArrayList<Integer> getAssigned( Collection<Integer> children,Set<Integer> keys){
		ArrayList<Integer> assigned = new ArrayList<Integer>();
		for(int key : keys)
			if(children.contains(key))
				assigned.add(key);
		return assigned;
	} 
	private HashSet<Integer> getPath(HashSet<Integer>children, Set<Integer> keys, int target){
		HashSet<Integer> parents = new  HashSet<Integer>();
		
		for(int child : children ){
    		if(keys.contains(child)){
    			positionsToKeep.add(child);
    		}
    	int id = Integer.parseInt(ph.getNode(String.valueOf(child)).getParent().getName());
    	if(id != target){
    		parents.add(id);
    		}
		}
		return  parents;
}
	public ArrayList<Integer> getAllStrains(int target, Set<Integer> keys){
		HashSet<Integer> children = new HashSet<Integer>();
		HashSet<Integer> parents = new HashSet<Integer>();
		int maxDepth = 0;
	    for(PhylogenyNode test : ph.getNode(String.valueOf(target)).getAllExternalDescendants()){
	      if(maxDepth < test.calculateDepth()){
	    	  maxDepth = test.calculateDepth();
	    	  }
	      	 children.add(Integer.parseInt(test.getName()));
	      	 int id = Integer.parseInt(test.getParent().getName());
	      	 if(id!=target)
	      		 parents.add(id);
	      	 
	    }	
	    positionsToKeep.addAll(getAssigned(children,keys));
	    for(int i = 0;i< (maxDepth-ph.getNode(String.valueOf(target)).calculateDepth());i++){
	    	parents = getPath(parents, keys,target);
	    	if(parents.size() == 0)
	    		break;
	    }
	   ArrayList<Integer> results = new ArrayList<Integer>();
			   results.addAll(positionsToKeep);
	   
		return results;
	}
		public ArrayList<Integer>  getParents(int target){
			ArrayList<Integer> ids = new ArrayList<Integer>();
			ids.add(target);
			int id = target;
			for(int i = 0; i<ph.getNode(String.valueOf(target)).calculateDepth();i++){
				PhylogenyNode t = ph.getNode(String.valueOf(id)).getParent();
				id = Integer.parseInt(t.getName());
				ids.add(id);
			}
			
			return ids;
		}
}
