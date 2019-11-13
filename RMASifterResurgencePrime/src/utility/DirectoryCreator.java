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
		
		switch(behave) {
			case CRAWL:
				if(wantMeganSummaries){
					new File(outDir+"/MeganSummaries/").mkdirs();
				}
				new File(outDir+"/crawlResults/").mkdirs();
				new File(outDir+"/crawlResults/damageMismatch/").mkdirs();
				new File(outDir+"/crawlResults/readDist/").mkdirs();
				new File(outDir+"/crawlResults/editDistance/").mkdirs();
				new File(outDir+"/crawlResults/percentIdentity/").mkdirs();
				new File(outDir+"/crawlResults/reads/").mkdirs();
				new File(outDir+"/crawlResults/coverage/").mkdirs();
			break;
			
			case NON:
				new File(outDir+"/default/").mkdirs();
				new File(outDir+"/default/"+"/readDist/").mkdirs(); //TODO could break potentially on Windows systems
				new File(outDir+"/default/"+"/coverage/").mkdirs();
				new File(outDir+"/default/"+"/editDistance/").mkdirs();
				new File(outDir+"/default/"+"/percentIdentity/").mkdirs();
				new File(outDir+"/default/"+"/damageMismatch/").mkdirs();
				new File(outDir+"/default/"+"/filterInformation/").mkdirs();
				if(hits)
					new File(outDir+"/default/"+"/alignments/").mkdirs();
				if(reads)
					new File(outDir+"/default/"+"/reads/").mkdirs();
			break;	
			
			case ANCIENT:
				new File(outDir+"/ancient/").mkdirs();
				new File(outDir+"/ancient/"+"/readDist/").mkdirs(); 
				new File(outDir+"/ancient/"+"/editDistance/").mkdirs();
				new File(outDir+"/ancient/"+"/percentIdentity/").mkdirs();
				new File(outDir+"/ancient/"+"/coverage/").mkdirs();
				new File(outDir+"/ancient/"+"/damageMismatch/").mkdirs();
				new File(outDir+"/ancient/"+"/filterInformation/").mkdirs();
				if(hits)
					new File(outDir+"/ancient/"+"/alignments/").mkdirs();
				if(reads)
					new File(outDir+"/ancient/"+"/reads/").mkdirs();
			break;	
			
			case NON_ANCIENT:
				new File(outDir+"/default/").mkdirs();
				new File(outDir+"/ancient/"+"/readDist/").mkdirs(); 
				new File(outDir+"/ancient/"+"/editDistance/").mkdirs();
				new File(outDir+"/ancient/"+"/percentIdentity/").mkdirs();
				new File(outDir+"/ancient/"+"/damageMismatch/").mkdirs();
				new File(outDir+"/ancient/"+"/filterInformation/").mkdirs();
				new File(outDir+"/ancient/"+"/coverage/").mkdirs();
				
				new File(outDir+"/default/"+"/readDist/").mkdirs(); 
				new File(outDir+"/default/"+"/editDistance/").mkdirs();
				new File(outDir+"/default/"+"/percentIdentity/").mkdirs();
				new File(outDir+"/default/"+"/damageMismatch/").mkdirs();
				new File(outDir+"/default/"+"/filterInformation/").mkdirs();
				new File(outDir+"/default/"+"/coverage/").mkdirs();
				if(hits){
					new File(outDir+"/default/"+"/alignments/").mkdirs();
					new File(outDir+"/ancient/"+"/alignments/").mkdirs();
				}
				if(reads){
					new File(outDir+"/default/"+"/reads/").mkdirs();
					new File(outDir+"/ancient/"+"/reads/").mkdirs();
				}
				break;
			
			case SRNA: 
				new File(outDir+"/default/").mkdirs();
				new File(outDir+"/ancient/"+"/readDist/").mkdirs(); 
				new File(outDir+"/ancient/"+"/editDistance/").mkdirs();
				new File(outDir+"/ancient/"+"/percentIdentity/").mkdirs();
				new File(outDir+"/ancient/"+"/damageMismatch/").mkdirs();
				new File(outDir+"/ancient/"+"/filterInformation/").mkdirs();
				new File(outDir+"/ancient/"+"/coverage/").mkdirs();
				
				new File(outDir+"/default/"+"/readDist/").mkdirs(); 
				new File(outDir+"/default/"+"/editDistance/").mkdirs();
				new File(outDir+"/default/"+"/percentIdentity/").mkdirs();
				new File(outDir+"/default/"+"/damageMismatch/").mkdirs();
				new File(outDir+"/default/"+"/filterInformation/").mkdirs();
				new File(outDir+"/default/"+"/coverage/").mkdirs();
				if(hits){
					new File(outDir+"/default/"+"/alignments/").mkdirs();
					new File(outDir+"/ancient/"+"/alignments/").mkdirs();
				}
				if(reads){
					new File(outDir+"/default/"+"/reads/").mkdirs();
					new File(outDir+"/ancient/"+"/reads/").mkdirs();
				}
				break;
				case ASSIGNMENT:{
					
					break;
				}
				default:
					System.err.println("Filter no longer supprted use parameter -h for accepted values");
				break;	
			}
		}

}
