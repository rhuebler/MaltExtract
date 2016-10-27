package utility;

import java.math.BigDecimal;
import java.math.RoundingMode;
import java.util.ArrayList;

import RMAAlignment.Alignment;

public class removedStaitsitcs {
	private static BigDecimal getMean(ArrayList<Integer> data)
	{
	    float sum = 0;
	    for(int a : data)
	        sum += a;
	    return new BigDecimal(sum/data.size());
	}
	private static BigDecimal getVariance(ArrayList<Integer> data, BigDecimal mean){;
		BigDecimal temp = new BigDecimal(0);
	    for(double d : data){
	    	BigDecimal a=new BigDecimal(d);
	    	temp = temp.add((a.subtract(mean).multiply(a.subtract(mean))));		
	    }
	    return temp.divide(new BigDecimal(data.size()-1), RoundingMode.HALF_EVEN);
	    }

	private static BigDecimal getStdDev(BigDecimal variance)
	{
	    return new BigDecimal(Math.sqrt(variance.doubleValue()));
	}

	private static BigDecimal getMedian(ArrayList<Integer> data) 
	{	 double pos1 = Math.floor((data.size() - 1.0) / 2.0);
	     double pos2 = Math.ceil((data.size() - 1.0) / 2.0);
	     if (pos1 == pos2 ) {
	        return new BigDecimal(data.get((int)pos1));
	     } else {
	        return new BigDecimal(((data.get((int)pos1) + data.get((int)pos2)) / 2.0));
	     }
	}
//	length = cEnd - cStart + 1;
//	
//	
//	if(i==0){
//		int nStart = 0;
//		Alignment next = input.get(i+1);
//		if(next.isReversed()){
//			nStart = next.getEnd();	
//		}else{
//			nStart = next.getStart();
//		}
//		
//		int a = (cEnd-nStart);
//		if(a > 0)
//			soli = length - a;
//		else if(a == 0)
//			soli = length - 1;
//		else 
//			soli = length;
//		distance.add(nStart-cStart);
//		
//	}else if(i < input.size() - 1){
//		int nStart = 0;
//		int pEnd = 0;
//		Alignment next = input.get(i+1);
//		Alignment previous = input.get(i-1);
//		
//		if(next.isReversed()){
//			nStart = next.getEnd();	
//		}else{
//			nStart = next.getStart();		
//		}
//		
//		if(previous.isReversed()){
//			pEnd = previous.getStart();		
//		}else{
//			pEnd = previous.getEnd();		
//		}
//		int a = cEnd - nStart; 
//		int b =	pEnd - cStart;
//		if(a<0 && b<0)
//			soli = length;
//		else if(a<0 && b>0)
//			soli = length - a;
//		else if(a>0 && b<0)
//			soli = length - b;
//		if(a==0)soli -= 1;// test if either distance is zero and one position less can be considered unique 
//		if(b==0)soli -= 1;
//		distance.add(nStart-cStart);
//	
//	}else{
//		Alignment previous = input.get(i-1);
//		int pEnd = 0;
//		
//		if(previous.isReversed()){
//			pEnd = previous.getStart();		
//		}else{
//			pEnd = previous.getEnd();		
//		}
//		int a = pEnd - cStart;
//		if(a>0)
//			soli = length-a;
//			else if(a==0)
//				soli = length-1;
//			else if(a<0)
//				soli= length;
//	}	
//	if(soli <= 0 || soli >= length)
//		soli = length;
//	unique += soli;
//	possible += length;
//	i++;
//}
//distance.sort(null);
}
