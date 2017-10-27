package Read;

import RMAAlignment.Alignment;

public class Read {
	private String readName;
	private int readLength;
	private Alignment top;
	private String readSequence;
	public Read(String name, int length, String sequence){
		readName = name;
		readLength = length;
		readSequence = sequence;
	}
	public void setTopAlignment(Alignment alignment){
		this.top = alignment;
		}
	 public int getReadLength(){
		 return readLength;
	 }
	public String getReadName(){
		 return  readName;
	 }
	public Alignment getTopAlignment(){
		 return this.top;
	}
	public String getReadSequence(){
		return readSequence;
	}
}
