package RMAAlignment;

import java.math.BigDecimal;
import java.math.RoundingMode;
import java.util.ArrayList;
import java.util.List;
/**
 * This class is used to compute some Statistics from a List of Alignments. Alignments should come from CompositionMap and must have their duplicates marked
 * @author huebler
 *
 */
public class AlignmentStatistics {
	private ArrayList<Alignment> currentList;
	public AlignmentStatistics(ArrayList<Alignment> list){
		this.currentList = list;
	}
	//TODO add ALexanders moving window function 
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
	
	private ArrayList<Alignment> removeDuplicates(ArrayList<Alignment> input){
		ArrayList<Alignment> positionsToKeep = new ArrayList<Alignment>();
		for(Alignment al : input){
			if(!al.isDuplicate())
				positionsToKeep.add(al);
		}
			return positionsToKeep;
		
	}
	
	// process best list of start positions
	public List<Double> getStatistics(){
		ArrayList<Alignment> input = removeDuplicates(this.currentList);
		if(input.size() >= 2){
		ArrayList<Integer> distance = new ArrayList<Integer>();
		List<Double> results = new ArrayList<Double>();
		
		int i = 0;//TODO maybe there is a better solution also is there an logical error in here 
		double unique = 0;
		double possible = 0;
		while(i<input.size()){
			//calculate unique positions per read plus average distance between current and next read
			Alignment current = input.get(i);
			int length = 0;
			int cStart = 0;
			int cEnd = 0;
			int soli = 0;
			
			if(current.isReversed()){
				cStart = current.getEnd();
				cEnd = current.getStart();		
			}else{
				cStart = current.getStart();
				cEnd = current.getEnd();		
			}
			
			length = cEnd - cStart + 1;
			
			
			if(i==0){
				int nStart = 0;
				Alignment next = input.get(i+1);
				if(next.isReversed()){
					nStart = next.getEnd();	
				}else{
					nStart = next.getStart();
				}
				
				int a = (cEnd-nStart);
				if(a > 0)
					soli = length - a;
				else if(a == 0)
					soli = length - 1;
				else 
					soli = length;
				distance.add(nStart-cStart);
				
			}else if(i < input.size() - 1){
				int nStart = 0;
				int pEnd = 0;
				Alignment next = input.get(i+1);
				Alignment previous = input.get(i-1);
				
				if(next.isReversed()){
					nStart = next.getEnd();	
				}else{
					nStart = next.getStart();		
				}
				
				if(previous.isReversed()){
					pEnd = previous.getStart();		
				}else{
					pEnd = previous.getEnd();		
				}
				int a = cEnd - nStart; 
				int b =	pEnd - cStart;
				if(a<0 && b<0)
					soli = length;
				else if(a<0 && b>0)
					soli = length - a;
				else if(a>0 && b<0)
					soli = length - b;
				if(a==0)soli -= 1;// test if either distance is zero and one position less can be considered unique 
				if(b==0)soli -= 1;
				distance.add(nStart-cStart);
			
			}else{
				Alignment previous = input.get(i-1);
				int pEnd = 0;
				
				if(previous.isReversed()){
					pEnd = previous.getStart();		
				}else{
					pEnd = previous.getEnd();		
				}
				int a = pEnd - cStart;
				if(a>0)
					soli = length-a;
					else if(a==0)
						soli = length-1;
					else if(a<0)
						soli= length;
			}
				
			if(soli <= 0 || soli >= length)
				soli = length;
			unique += soli;
			possible += length;
			i++;
		}
		distance.sort(null);
		BigDecimal mean = getMean(distance);
		BigDecimal median = getMedian(distance);
		BigDecimal variance = getVariance(distance,mean);
		BigDecimal std = getStdDev(variance);
		results.add(mean.doubleValue());
		results.add(median.doubleValue());
		results.add(variance.doubleValue());
		results.add(std.doubleValue());
		results.add(unique/(possible)); //TODO consider taking actual average match length maybe calculate numbers without removing stacked reads 
		results.add((double)input.size());
		results.add((double) this.currentList.size());
		results.add((double)input.get(0).getReferenceLength());
		return results;
		}else{
		List<Double> results =	new ArrayList<Double>();
		for(int i = 0; i<8;  i++)
			results.add(0.0);
		return results;
		}
	}
}
