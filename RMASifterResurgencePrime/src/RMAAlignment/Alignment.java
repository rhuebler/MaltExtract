package RMAAlignment;

import java.util.Comparator;
/**
 * This class Represents an MatchBlock from Megan but has more Slots available
 * @author huebler
 *
 */
// readnamen, Node, wieviele Reads fallen da hon wo sie hinfallen und wieviele von denen target reads sind echt  
public class Alignment {
	private String readName;
	private double pIdent;
	private String query;
	private String reference;
	private String alignment;
	private String referenceName;
	private String accessionNumber;
	private int start;
	private int end;
	private int mLength;
	private int referenceLength;
	private int numGaps;
	private String strand;
	private boolean fivePrimeDamage;
	private boolean reversed = false;
	private boolean duplicate = false;
	private int editDistance;
	private int readLength;
	private double score;
	private String sequence;
	private String text;
	private boolean stacked = false;
	// getter
public boolean isStacked(){
	return this.stacked;
}	
public String getText(){
	return this.text;
}
public String getSequence(){
	return this.sequence;
}
public String getAccessionNumber(){
	return this.accessionNumber;
}	
public double getScore(){
	return this.score;
}	
public int getReadLength(){
	return this.readLength;
}	
public int getEditDistance(){
	return this.editDistance;
}
public String getReadName(){
	return this.readName;
}	
public double getPIdent(){
	return this.pIdent;
}
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
 public void setStacked(boolean b){
	 this.stacked = b;
 }
 public void setText(String s){
	 this.text = s;
 }
 public void setSequence(String s){
	 this.sequence = s;
 }
 public void setAcessionNumber(String s){
	 this.accessionNumber =s;
 }
public void setReadLength(int d){
	this.readLength = d;
}
public void setReadName(String s){
		this.readName = s;
	}	
 public void setPIdent(double d){
		this.pIdent = d;
	}
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
 
 private boolean misMatchFivePrime(int i){
	 if((reference.charAt(i) == 'C' && query.charAt(i) == 'T')||
		(reference.charAt(i) == 'c' && query.charAt(i) == 't')){
		 		return true;
	 }else {
	 		return false;
	 }
 }
 private boolean misMatchThreePrime(int i){
	 if((reference.charAt(i) == 'G' && query.charAt(i) == 'A')||
		(reference.charAt(i) == 'g' && query.charAt(i) == 'a')){
			 return true; 		
	 }else{
		 return false;
	 }
 }
 private void calculateEditDistance(String sequence, String reference){
	 int len1 = sequence.length();
		int len2 = reference.length();
		 int[][] d = new int[len1 + 1][len2 + 1];
		 for(int i = 0; i <= len1; i++)
			 d[i][0] = i;
		 for(int j = 0; j <= len2; j++)
			 d[0][j] = j;
		 for(int i = 0; i < len1; i++)
			 for(int j = 0; j < len2; j++){
				 if(sequence.charAt(i) == reference.charAt(j)){
					d[i+1][j+1] = d[i][j]; 
				 }else{
					int replace  = d[i][j]+1;
					int insert = d[i + 1][j]+1;
					int delete = d[i][j + 1] + 1;
					d[i+1][j+1] = Math.min(replace, 
									 Math.min(insert,
									 delete));
				 }
			 }	 
	 this.editDistance = d[sequence.length()][reference.length()];
 }
 
	 public void processText(){
		 String alignment="";
		 String query ="";
		 String reference = "";
		 boolean first = true;
		 for(String line : text.split("\n")){
			 if(line.startsWith(">"))
				 setReferenceName(line.trim());
			 if(line.contains("Query")){
				 for(String frag:line.split(":|\\s")){
					 if(frag.trim().matches("[ATGCNatgcn-]+")){
						 query+=(frag.trim());
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
			 }	 
			 if(line.contains("Sbjct")){
				
				 for(String frag:line.split(":|\\s")){
					 if(frag.trim().matches("[ATGCNatgcn-]+")){
						 reference += (frag.trim());	 
					 }
					 if(frag.trim().matches("\\d+")){
						 if(first){
							setStart(Integer.parseInt(frag.trim()));
							first = false;
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
		 setQuery(query);
		 setAlignment(alignment);
		 setReference(reference);
		 if(start > end){
			 setReversed(true);
		 }
		 for(int k = 0; k < 5; k++ ){ //test first and last 5 positions for g->c mismatch
			 if(reference.length() > 1+k*2 && query.length() > 1+k*2 
					 && reference.length() == query.length()){
				 if(misMatchFivePrime(k) || misMatchThreePrime(query.length()-k-1)){ // set damage true if there is a C->T sub at eithter end
					setFivePrimeDamage(true);
					break;// break out of loop
				}
		 	}
		 }
		 calculateEditDistance(query, reference);
	 }

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