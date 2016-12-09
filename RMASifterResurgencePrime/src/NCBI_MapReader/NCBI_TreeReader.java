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
		private ArrayList<Integer> getStrains(ArrayList<Integer>children, Set<Integer> keys){
			ArrayList<Integer> first = new  ArrayList<Integer>();
			ArrayList<Integer> second = new  ArrayList<Integer>();
	    for(int child : children ){
	    		for(PhylogenyNode grandChild : ph.getNode(String.valueOf(child)).getDescendants()){
	    			first.add(Integer.parseInt(grandChild.getName()));
	    			for(PhylogenyNode greatGrandChild:grandChild.getDescendants()){
	    			second.add(Integer.parseInt(greatGrandChild.getName()));
	    			}
	    		}	
	    		
	    }	
	    positionsToKeep.addAll(getAssigned(children,keys));
	    first.addAll(second);
	    return  getAssigned(first,keys);
	    
	}
		public ArrayList<Integer> getAllStrains(int target, Set<Integer> keys){
			ArrayList<Integer> children = new ArrayList<Integer>();
			int maxDepth = 0;
		    for(PhylogenyNode test : ph.getNode(String.valueOf(target)).getAllExternalDescendants()){
		      if(maxDepth < test.calculateDepth())
		    	  maxDepth = test.calculateDepth();
		    }
		    for(PhylogenyNode test : ph.getNode(String.valueOf(target)).getDescendants()){
		    	children.add(Integer.parseInt(test.getName()));
		    }
		    
		    ArrayList<Integer> positions = new  ArrayList<Integer>();
		    positions.addAll(children);
		    for(int i = 0;i< (maxDepth-ph.getNode(String.valueOf(target)).calculateDepth());i++){
		    	positions = getStrains(positions, keys);
		    	if(positions.size() == 0)
		    		break;
		    }
		    positionsToKeep.addAll(getAssigned(positions,keys));
		    //System.out.println(positionsToKeep.size());
		   children.clear();
		    children.addAll(positionsToKeep);
			return children;
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
