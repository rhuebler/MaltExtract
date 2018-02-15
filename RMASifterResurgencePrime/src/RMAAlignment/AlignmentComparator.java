package RMAAlignment;
/**
 * Comparator for RMA6 Alignmet objects necessary for sorting to allow PCR duplicate removal
 * @author huebler
 */
import java.util.Comparator;

public class AlignmentComparator implements Comparator<Alignment> // comparator for alignment class that takes into consideration if the read
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