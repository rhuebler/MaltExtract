package utility;

import java.io.File;

import behaviour.Filter;

public class DirectoryCreator {
/**
 * @author huebler
 * Create all directories associated with RMAextractor depending on input parameters
 * @param behave
 * @param outDir
 * @param hits
 * @param crawl
 * @param reads
 */
	public void process(Filter behave, String outDir, boolean hits,boolean reads, boolean wantMeganSummaries) {
		if(behave != Filter.CRAWL && wantMeganSummaries){
			new File(outDir+"/MeganSummaries/").mkdirs();
		}
		if(behave == Filter.CRAWL){
			new File(outDir+"/crawlResults/").mkdirs();
			new File(outDir+"/crawlResults/damageMismatch/").mkdirs();
			new File(outDir+"/crawlResults/readDist/").mkdirs();
			new File(outDir+"/crawlResults/editDistance/").mkdirs();
			new File(outDir+"/crawlResults/percentIdentity/").mkdirs();
			new File(outDir+"/crawlResults/reads/").mkdirs();
		}
		if(behave == Filter.NON){
			new File(outDir+"/default/").mkdirs();
			new File(outDir+"/default/"+"/readDist/").mkdirs(); //TODO could break potentially on Windows systems
			new File(outDir+"/default/"+"/editDistance/").mkdirs();
			new File(outDir+"/default/"+"/percentIdentity/").mkdirs();
			new File(outDir+"/default/"+"/damageMismatch/").mkdirs();
			if(hits)
				new File(outDir+"/default/"+"/alignments/").mkdirs();
			if(reads)
				new File(outDir+"/default/"+"/reads/").mkdirs();
		}else if(behave == Filter.ANCIENT){
			new File(outDir+"/ancient/").mkdirs();
			new File(outDir+"/ancient/"+"/readDist/").mkdirs(); 
			new File(outDir+"/ancient/"+"/editDistance/").mkdirs();
			new File(outDir+"/ancient/"+"/percentIdentity/").mkdirs();
			new File(outDir+"/ancient/"+"/damageMismatch/").mkdirs();
			if(hits)
				new File(outDir+"/ancient/"+"/alignments/").mkdirs();
			if(reads)
				new File(outDir+"/ancient/"+"/reads/").mkdirs();
		}else if(behave == Filter.NONDUPLICATES){
			new File(outDir+"/nonDuplicates/").mkdirs();
			new File(outDir+"/nonDuplicates/"+"/readDist/").mkdirs(); 
			new File(outDir+"/nonDuplicates/"+"/editDistance/").mkdirs();
			new File(outDir+"/nonDuplicates/"+"/percentIdentity/").mkdirs();
			new File(outDir+"/nonDuplicates/"+"/damageMismatch/").mkdirs();
			if(hits)
				new File(outDir+"/nonDuplicates/"+"/alignments/").mkdirs();
			if(reads)
				new File(outDir+"/nonDuplicates/"+"/reads/").mkdirs();
		}else if(behave == Filter.ALL){
			new File(outDir+"/ancientNonDuplicates/").mkdirs();
			new File(outDir+"/ancientNonDuplicates/"+"/readDist/").mkdirs(); 
			new File(outDir+"/ancientNonDuplicates/"+"/editDistance/").mkdirs();
			new File(outDir+"/ancientNonDuplicates/"+"/percentIdentity/").mkdirs();
			new File(outDir+"/ancientNonDuplicates/"+"/damageMismatch/").mkdirs();
			if(hits)
				new File(outDir+"/ancientNonDuplicates/"+"/alignments/").mkdirs();
			if(reads)
				new File(outDir+"/ancientNonDuplicates/"+"/reads/").mkdirs();
		}else if(behave == Filter.NON_ANCIENT){
			new File(outDir+"/default/").mkdirs();
			new File(outDir+"/ancient/"+"/readDist/").mkdirs(); 
			new File(outDir+"/ancient/"+"/editDistance/").mkdirs();
			new File(outDir+"/ancient/"+"/percentIdentity/").mkdirs();
			new File(outDir+"/ancient/"+"/damageMismatch/").mkdirs();
			new File(outDir+"/ancient/"+"/FilterInformation/").mkdirs();
			
			new File(outDir+"/default/"+"/readDist/").mkdirs(); 
			new File(outDir+"/default/"+"/editDistance/").mkdirs();
			new File(outDir+"/default/"+"/percentIdentity/").mkdirs();
			new File(outDir+"/default/"+"/damageMismatch/").mkdirs();
			new File(outDir+"/default/"+"/FilterInformation/").mkdirs();
			if(hits){
				new File(outDir+"/default/"+"/alignments/").mkdirs();
				new File(outDir+"/ancient/"+"/alignments/").mkdirs();
			}
			if(reads){
				new File(outDir+"/default/"+"/reads/").mkdirs();
				new File(outDir+"/ancient/"+"/reads/").mkdirs();
			}
			}
	}

}
