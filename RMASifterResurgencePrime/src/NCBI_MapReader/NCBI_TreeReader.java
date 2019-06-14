package NCBI_MapReader;
import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.net.URL;
import java.net.URLConnection;
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
 * Unfortunately Malt does only store taxIds of nodes that have at least one assigned 
 * read therefore the number of nodes cannot be reduced at every level for the risk of losing 
 * a branch. Instead a more crawling like approach is currently needed that allows to skip one taxonomic Level or two 
 * @author huebler
 *
 */

public class NCBI_TreeReader {
	/**
	 * @param String directory/to/file
	 * @throws none thrown all caught
	 * @return return phylogenetic tree object
	 */
	
	private String treName= "";
	int target;
	private HashSet<Integer> positionsToKeep = new HashSet<Integer>();
	private Phylogeny ph;
	public NCBI_TreeReader(){
		processFromWeb();
	}
	
	public void processFromWeb() {
		try{
			
			String location ="https://raw.githubusercontent.com/rhuebler/HOPS/SHH/AMPS/Resources/ncbi.tre";
			URLConnection conn = new URL(location).openConnection();
			 conn.setConnectTimeout(90*1000);
			 conn.setReadTimeout(90*1000);
			   try (InputStream in = conn.getInputStream()) {
				   InputStreamReader reader = new  InputStreamReader(in);
				   BufferedReader buffered = new BufferedReader(reader);
				   String line;
				   while((line = buffered.readLine())!=null) {
					
					   this.ph = Phylogeny.createInstanceFromNhxString(line);
			       }
				 buffered.close();
			   }catch(Exception e) {
				   e.printStackTrace();
			    }
		}catch(IOException io) {
			io.printStackTrace();
		}
	}
	
	public NCBI_TreeReader(String path){
		if(!path.endsWith("/"))
			path+="/";
		this.treName = path + "ncbi.tre";// if path is provided use path
		try{
			Scanner in = new Scanner(new File(treName));
			String line = in.nextLine();
			in.close();
		    this.ph = Phylogeny.createInstanceFromNhxString(line);
		    }catch(IOException io){
		    	io.printStackTrace();
		    }
	}
	// get files by downloading 
	public NCBI_TreeReader(NCBI_TreeReader copyInstance){
		this.ph = copyInstance.getPhylogeny();
	}
	
	private Phylogeny getPhylogeny(){
		return ph;
	}
	
	private ArrayList<Integer> getAssigned( Collection<Integer> children,Set<Integer> keys){// get assinged nodes in files
		ArrayList<Integer> assigned = new ArrayList<Integer>();
		if(keys==null)
			System.err.println("Danger empty keys in File");
		for(int key : keys)
			if(children.contains(key))
				assigned.add(key);
		return assigned;
	} 
	
	public ArrayList<Integer> getAllStrains(int target, Set<Integer> keys){// get all strains of a node
		positionsToKeep.clear();
		PhylogenyNode query = ph.getNode(String.valueOf(target));
	    for(PhylogenyNode leaf : query.getAllExternalDescendants()){
	    	HashSet<Integer> children = new HashSet<Integer>();
	    	children.add(Integer.parseInt(leaf.getName()));
	      for(int i = 0;i< (leaf.calculateDepth()-query.calculateDepth());i++){
	    	leaf = leaf.getParent();
	    	int id = Integer.parseInt(leaf.getName());
	    	if(id != target){
	    		children.add(id);
	    	}
	      }
	      positionsToKeep.addAll(getAssigned(children,keys));
	    }	
	    ArrayList<Integer> results = new ArrayList<Integer>();
		results.addAll(positionsToKeep);
		return results;
	}   
		public  ArrayList<Integer> getLeavesIds(){
			positionsToKeep.clear();
			PhylogenyNode query = ph.getNode(0);
			for(PhylogenyNode leaf : query.getAllExternalDescendants()){
				positionsToKeep.add(Integer.parseInt(leaf.getName()));
			}
			ArrayList<Integer> results = new ArrayList<Integer>();
			   results.addAll(positionsToKeep);
			   return results;
		}
	    public ArrayList<Integer> getLowestContainedIds(Set<Integer> keys){// get all strains of a node
	    	positionsToKeep.clear();
			PhylogenyNode query = ph.getNode(0);
		    for(PhylogenyNode leaf : query.getAllExternalDescendants()){
		    	if(keys.contains(Integer.parseInt(leaf.getName()))) {
		    		positionsToKeep.add(Integer.parseInt(leaf.getName()));
		    	}else {
		    		 for(int i = 0;i< (leaf.calculateDepth()-query.calculateDepth());i++){
		    		    	leaf = leaf.getParent();
		    		    	int id = Integer.parseInt(leaf.getName());
		    		    	if(keys.contains(id)){
		    		    		positionsToKeep.add(id);
		    		    		break;
		    		    	}
		    		      }
		    	}	

		     
		    }	
	   ArrayList<Integer> results = new ArrayList<Integer>();
			   results.addAll(positionsToKeep);
	   
		return results;
	   
	}
	    
		public ArrayList<Integer>  getParents(int target){// get parents of leaves
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
