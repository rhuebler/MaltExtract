package RMAAlignment;

import java.util.Comparator;

public class Alignment {
	private String query;
	private String reference;
	private String alignment;
	private String referenceName;
	private int start;
	private int end;
	private int mLength;
	private int referenceLength;
	private int numGaps;
	private String strand;
	private boolean fivePrimeDamage;
	private boolean reversed = false;
	private boolean duplicate = false;
	// getter
public boolean isReversed(){
	return this.reversed;
}	
public boolean isDuplicate(){
	return this.duplicate;
}	
public String getReferenceName(){
	return this.referenceName;}	

public int getNumGaps(){	
	return this.numGaps;}	
	
public boolean getFivePrimeDamage(){
	return this.fivePrimeDamage;}

public String getStrand(){
	return this.strand;}

public String getQuery(){
	 return this.query;}

 public String getReference(){
	 return this.reference;}
 
 public String getAlignment(){
	 return this.alignment;}
 
 public int getStart(){
	 return this.start;}
 
 public int getEnd(){
	 return this.end;}
 
 public int getMlength(){
	 return this.mLength;}
 
 public int getReferenceLength(){
	 return this.referenceLength;}
 
 // setters
 public void setDuplicate(boolean b){ //works
	 this.duplicate = b;
 }
 private void setReversed(boolean b){
	 this.reversed = b;
 }
 private void setReferenceName(String s){
	 this.referenceName = s;}
 
 private void setQuery(String query){
	 this.query = query;}
 
 private void setReference(String reference){
	 this.reference = reference;}
 
 private void setAlignment(String alignment){
	 this.alignment = alignment;}
 
 private void setStart(int start){
	 this.start = start;}
 
 private void setEnd(int end){
	 this.end = end;}
 
 private void setMLength(int length){
	 this.mLength = length;}
 
 private void setStrand(String s) {
		this.strand = s;}
 
private void setFivePrimeDamage(boolean b) {
		this.fivePrimeDamage = b;}

private void setReferenceLength(int k) {
	this.referenceLength=k;}

 private void setNumGaps(int k){
	 this.numGaps = k;}
 
 private boolean ctMisMatch(int i){
	 if((this.reference.charAt(i) == 'C' && this.query.charAt(i) == 'T')
				|| this.reference.charAt(i) == 'c' && this.query.charAt(i) == 't')
		 		return true;
	 else return false;
 }
 
 // process matchblock and save information in Alignment object
 public void processText(String[] text){
	 String alignment="";
	 for(String line : text){
		 if(line.startsWith(">"))
			 setReferenceName(line.trim());
		 if(line.contains("Query")){
			 for(String frag:line.split("\\s")){
				 if(frag.trim().matches("[ATGCNatgcn-]+")){
					 setQuery(frag.trim());
					 setMLength(frag.trim().length());
					 }
			 }
		 }
		 if(line.matches("[\\s\\|]+")){
			 
			 for(char c : line.substring(line.length()-mLength,line.length()).toCharArray())
			 { if(c=='|'){
					 alignment += 1;
				 }else{
					 alignment += 0;
				 }
		    }
			setAlignment(alignment); 
		 }	 
		 if(line.contains("Sbjct")){
			 int i =0;
			 for(String frag:line.split("\\s")){
				 if(frag.trim().matches("[ATGCatgc]+")){
					 setReference(frag.trim());
					if(ctMisMatch(0) || ctMisMatch(this.query.length()-1)) // set damage true if there is a C->T sub at eithter end
						setFivePrimeDamage(true);
				 	}
				 if(frag.trim().matches("\\d+")){
					 if(i==0){
						setStart(Integer.parseInt(frag.trim()));
						 i++;
					 }else{
						 setEnd(Integer.parseInt(frag.trim()));
					 }	 
				 }
			 }//for
	    }//outer if
		if(line.contains("Strand"))
			for(String frag : line.split("=")){
				if(frag.trim().matches("\\w+\\s/\\s\\w+"))
					setStrand(frag.trim());
			}
		if(line.contains("Length"))
			for(String frag : line.split("="))
				if(frag.trim().matches("\\d+"))
					setReferenceLength(Integer.parseInt(frag.trim()));
		if(line.contains("Gaps"))
			setNumGaps(Integer.parseInt((line.split(",")[1].split("=")[1].split("/")[0]).trim()));
	}// for
	 if(getStart() > getEnd()){
		 setReversed(true);
	 }
 }// method

}//classbody

class AlignmentComparator implements Comparator<Alignment> // comparator for alignment class that takes into consideration if the read
{															// if the Read is on the reverse Strand //TODO double check
  @Override public int compare( Alignment al1, Alignment al2 )
  {
	if(!al1.isReversed() && !al2.isReversed()){
			return al1.getStart() - al2.getStart();
		}else if(al1.isReversed() && !al2.isReversed()){
			return al1.getEnd() - al2.getStart();
		}else if(!al1.isReversed() && al2.isReversed()){
		 return al1.getStart() - al2.getEnd();
		 }else{
			 return al1.getEnd() - al2.getEnd();
		 }
	
  }
}