package RMAAlignment;

public class NOAOR {
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
