package NCBI_MapReader;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Scanner;
import java.util.Set;

import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;
/**
 * Process NCBI tre File with archeoptrix.jar
 * Process tree find child nodes to input taxID
 * only keep child IDs that are assigned in current RMA6 file 
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
	private ArrayList<Integer> getAssigned( ArrayList<Integer> children,Set<Integer> keys){
		ArrayList<Integer> assigned = new ArrayList<Integer>();
		for(int key : keys)
			if(children.contains(key))
				assigned.add(key);
		return assigned;
	} 
		private ArrayList<Integer> getStrains(ArrayList<Integer>children, int target, Set<Integer> keys){
			ArrayList<Integer> targets = new  ArrayList<Integer>();
	    for(int child : children ){ // works in principal but would have to analyze a lot of nodes for it to work.... is there a smaller tre file ?
	    		if(child != target)
	    		 targets.add(Integer.parseInt(ph.getNode(String.valueOf(child)).getParent().getName()));
	    }	
	    return  getAssigned(targets,keys);
	    
	}
		public ArrayList<Integer> getAllStrains(int target, Set<Integer> keys){// maybe it is more efficient to follow the root up
			ArrayList<Integer> children = new ArrayList<Integer>();
			int maxDepth = 0;
		    for(PhylogenyNode test : ph.getNode(String.valueOf(target)).getAllDescendants()){ // works in principal but would have to analyze a lot of nodes for it to work.... is there a smaller tre file ?
		      if(maxDepth < test.calculateDepth())
		    	  maxDepth = test.calculateDepth();
		      	children.add(Integer.parseInt(test.getName()));
		    }
		    children = getAssigned(children,keys);
		    ArrayList<Integer> positions = new  ArrayList<Integer>();
		    positions.addAll(children);
		    for(int i = 0;i< maxDepth-ph.getNode(String.valueOf(target)).calculateDepth();i++){
		    	positions = getStrains(positions, target, keys);
		    	children.addAll(positions);
		    }
			return getAssigned(children,keys);
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
