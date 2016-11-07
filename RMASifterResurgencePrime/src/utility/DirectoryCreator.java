package utility;

import java.io.File;

import behaviour.Filter;

public class DirectoryCreator {

	public void process(Filter behave, String outDir, boolean hits, boolean crawl) {
		if(behave !=Filter.NON_ANCIENT){
		new File(outDir+"/readDist/").mkdirs(); //TODO could break potentially on Windows systems
		new File(outDir+"/editDistance/").mkdirs();
		new File(outDir+"/percentIdentity/").mkdirs();
		new File(outDir+"/damageMismatch/").mkdirs();
		if(hits)
			new File(outDir+"/reads/").mkdirs();
		if(crawl)
			new File(outDir+"/crawlResults/").mkdirs();
		}else if(behave==Filter.NON_ANCIENT){
		new File(outDir+"/default/").mkdirs();
		new File(outDir+"/ancient/").mkdirs();
		new File(outDir+"/ancient/"+"/readDist/").mkdirs(); //TODO could break potentially on Windows systems
		new File(outDir+"/ancient/"+"/editDistance/").mkdirs();
		new File(outDir+"/ancient/"+"/percentIdentity/").mkdirs();
		new File(outDir+"/ancient/"+"/damageMismatch/").mkdirs();
		new File(outDir+"/default/"+"/readDist/").mkdirs(); //TODO could break potentially on Windows systems
		new File(outDir+"/default/"+"/editDistance/").mkdirs();
		new File(outDir+"/default/"+"/percentIdentity/").mkdirs();
		new File(outDir+"/default/"+"/damageMismatch/").mkdirs();
		if(hits)
			new File(outDir+"/ancient/"+"/reads/").mkdirs();
		}

	}

}
