package utility;

import java.io.File;

import behaviour.Filter;

public class DirectoryCreator {

	public void process(Filter behave, String outDir, boolean hits, boolean crawl) {
		if(crawl){
			new File(outDir+"/crawlResults/").mkdirs();
		}
		if(behave == Filter.NON){
			new File(outDir+"/default/").mkdirs();
			new File(outDir+"/default/"+"/readDist/").mkdirs(); //TODO could break potentially on Windows systems
			new File(outDir+"/default/"+"/editDistance/").mkdirs();
			new File(outDir+"/default/"+"/percentIdentity/").mkdirs();
			new File(outDir+"/default/"+"/damageMismatch/").mkdirs();
			if(hits)
				new File(outDir+"/default/"+"/reads/").mkdirs();
		}else if(behave == Filter.ANCIENT){
			new File(outDir+"/ancient/").mkdirs();
			new File(outDir+"/ancient/"+"/readDist/").mkdirs(); 
			new File(outDir+"/ancient/"+"/editDistance/").mkdirs();
			new File(outDir+"/ancient/"+"/percentIdentity/").mkdirs();
			new File(outDir+"/ancient/"+"/damageMismatch/").mkdirs();
			if(hits)
				new File(outDir+"/ancient/"+"/reads/").mkdirs();
		}else if(behave == Filter.NONDUPLICATES){
			new File(outDir+"/nonDuplicates/").mkdirs();
			new File(outDir+"/nonDuplicates/"+"/readDist/").mkdirs(); 
			new File(outDir+"/nonDuplicates/"+"/editDistance/").mkdirs();
			new File(outDir+"/nonDuplicates/"+"/percentIdentity/").mkdirs();
			new File(outDir+"/nonDuplicates/"+"/damageMismatch/").mkdirs();
		}else if(behave == Filter.ALL){
			new File(outDir+"/ancientNonDuplicates/").mkdirs();
			new File(outDir+"/ancientNonDuplicates/"+"/readDist/").mkdirs(); 
			new File(outDir+"/ancientNonDuplicates/"+"/editDistance/").mkdirs();
			new File(outDir+"/ancientNonDuplicates/"+"/percentIdentity/").mkdirs();
			new File(outDir+"/ancientNonDuplicates/"+"/damageMismatch/").mkdirs();
			if(hits)
				new File(outDir+"/ancientNonDuplicates/"+"/reads/").mkdirs();
		}else if(behave == Filter.NON_ANCIENT){
			new File(outDir+"/default/").mkdirs();
			new File(outDir+"/ancient/"+"/readDist/").mkdirs(); 
			new File(outDir+"/ancient/"+"/editDistance/").mkdirs();
			new File(outDir+"/ancient/"+"/percentIdentity/").mkdirs();
			new File(outDir+"/ancient/"+"/damageMismatch/").mkdirs();
			
			new File(outDir+"/default/"+"/readDist/").mkdirs(); 
			new File(outDir+"/default/"+"/editDistance/").mkdirs();
			new File(outDir+"/default/"+"/percentIdentity/").mkdirs();
			new File(outDir+"/default/"+"/damageMismatch/").mkdirs();
			if(hits){
				new File(outDir+"/default/"+"/reads/").mkdirs();
				new File(outDir+"/ancient/"+"/reads/").mkdirs();
			}
			}
	}

}
