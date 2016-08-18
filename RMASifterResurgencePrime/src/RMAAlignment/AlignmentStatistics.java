package RMAAlignment;

import java.util.ArrayList;
import java.util.List;

public class AlignmentStatistics {
	private ArrayList<Alignment> currentList;
	public AlignmentStatistics(ArrayList<Alignment> list){
		this.currentList = list;
	}
	//TODO add ALexanders moving window function 
	private double getMean(ArrayList<Integer> data)
	{
	    double sum = 0.0;
	    for(int a : data)
	        sum += a;
	    return (sum/data.size());
	}
	private double getVariance(ArrayList<Integer> data, double mean){
	    double temp = 0;
	    for(double a :data)
	        temp += (mean-a)*(mean-a);
	    return temp/data.size();}

	private double getStdDev(double variance)
	{
	    return Math.sqrt(variance);
	}

	private double getMedian(ArrayList<Integer> data) 
	{	 double pos1 = Math.floor((data.size() - 1.0) / 2.0);
	     double pos2 = Math.ceil((data.size() - 1.0) / 2.0);
	     if (pos1 == pos2 ) {
	        return data.get((int)pos1);
	     } else {
	        return (data.get((int)pos1) + data.get((int)pos2)) / 2.0 ;
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
		ArrayList<Integer> distance = new ArrayList<Integer>();
		List<Double> results= new ArrayList<Double>();
		
		int i = 0;//TODO maybe there is a better solution also is there an logical error in here 
		double unique = 0;
		while(i<input.size()){
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
				
			}else if(i<input.size()-1){
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
			unique+=soli;
			i++;
		}
		
		double mean = getMean(distance);
		double median = getMedian(distance);
		double variance = getVariance(distance,mean);
		double std = getStdDev(variance);
		results.add(mean);
		results.add(median);
		results.add(variance);
		results.add(std);
		results.add(unique/(100.0*input.size())); //TODO consider taking actual average match length
		results.add((double)input.size());
		results.add((double) this.currentList.size());
		results.add((double)input.get(0).getReferenceLength());
		return results;
	}
}
