package RMAAlignment;

public class NOAOR {
	/**
	 * @author huebler
	 * This Class is used to sort the reference for the output of the the top ten references on a node
	 */
	private int size;
	private String reference;
	private int taxID;
	public NOAOR(int number, String name, int id){
		size = number;
		reference = name;	
		taxID = id;
	}
	public int getSize(){
		return size;
	}
	public String getReference(){
		return reference;
	}
	public int getTaxID(){
		return taxID;
	}
}
