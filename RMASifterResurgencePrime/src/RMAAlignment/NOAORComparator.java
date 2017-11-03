package RMAAlignment;

import java.util.Comparator;

public class NOAORComparator implements Comparator<NOAOR> // comparator for alignment class that takes into consideration if the read
{															// if the Read is on the reverse Strand //TODO double check
	  @Override public int compare(NOAOR n1, NOAOR n2 )
	  {
		  return n2.getSize()-n1.getSize();
	  }

}
