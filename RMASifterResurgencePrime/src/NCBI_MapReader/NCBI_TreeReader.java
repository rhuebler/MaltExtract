package NCBI_MapReader;
import java.io.File;
import java.util.ArrayList;
import java.util.Scanner;
import java.util.Set;
import java.io.IOException;    

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
	final static String treName= "/Users/huebler/Desktop/Project_RMASifter_The_RMAing/megan-ce-master/resources/files/ncbi.tre";// how to provide System resources
	int target;
	private Phylogeny ph;
	
	public NCBI_TreeReader() throws IOException{
		System.out.println("Setting up Phylogenetic Tree");
		Scanner in = new Scanner(new File(treName));
		String nhx = in.nextLine();
		in.close();
	    this.ph = Phylogeny.createInstanceFromNhxString(nhx);
	}
	private ArrayList<Integer> getAssigned( ArrayList<Integer> children,Set<Integer> keys){
		ArrayList<Integer> assigned = new ArrayList<Integer>();
		for(int key : keys)
			if(children.contains(key))
				assigned.add(key);
		return assigned;
	}
		public ArrayList<Integer> getStrains(int target, Set<Integer> keys){
	    ArrayList<Integer> children = new  ArrayList<Integer>();
	    for(PhylogenyNode test : ph.getNode(String.valueOf(target)).getDescendants()){ // works in principal but would have to analyze a lot of nodes for it to work.... is there a smaller tre file ?
	      
	    	children.add(Integer.parseInt(test.getName()));// so what if I follow Nodes to leaf and see if any thing is assinged there? 
	    	for (PhylogenyNode t : test.getDescendants()) {// should i solve this here? if in keys do something to get strain IDs and numbers 
	    		 children.add(Integer.parseInt(t.getName()));
			}
	    }	
	    return  getAssigned(children,keys);
	}
}
