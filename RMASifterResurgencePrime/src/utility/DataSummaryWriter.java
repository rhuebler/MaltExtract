package utility;

import java.io.File;
import java.io.IOException;
import java.nio.charset.Charset;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.logging.Level;
import java.util.logging.Logger;

import megan.core.DataTable;
import megan.core.MeganFile;
import megan.core.SampleAttributeTable;
import megan.core.SyncArchiveAndDataTable;
import megan.data.IConnector;

public class DataSummaryWriter {
	Logger warning;
	public DataSummaryWriter(Logger warning){
		this.warning = warning;
	}
	public void writeSummary(String directory,String fileName, String outDir) {
		// TODO Auto-generated method stub
		try {
			SampleAttributeTable sampleAttributeTable = new SampleAttributeTable();
			DataTable table = new DataTable();
			MeganFile file = new MeganFile();
			File f  = new File(directory+fileName);
			if(f.canRead() && f.canRead()){
				file.setFileFromExistingFile(f.getCanonicalPath(), true);
				file.setReadOnly(true);
				IConnector connector = file.getConnector();
				SyncArchiveAndDataTable.syncArchive2Summary(false, file.getFileName(), connector, table, sampleAttributeTable);
				if(!fileName.endsWith(".rma6"))
					fileName+=".rma6";
				ArrayList<String> output = new ArrayList<String>();
				output.add("@Creator\t"+ table.getCreator());
				output.add("@CreationDate\t"+ table.getCreationDate());
				output.add("@ContentType\t"+table.getContentType());
				output.add("@Names\t"+fileName);
				output.add("@BlastMode\t"+table.getBlastMode());
				output.add("@Uids\t"+table.getSampleUIds()[0]);
				output.add("@Sizes\t"+table.getSampleSizes()[0]);
				output.add("@TotalReads\t"+table.getTotalReads());
				output.add("@AdditionalReads\t"+table.getAdditionalReads());
				output.add("@Collapse\t"+"Taxonomy\t"+"-1\t");//table.getCollapsed("Taxonomy"));//what are you
				output.add("@Algorithm\t"+table.getAlgorithm("Taxonomy"));
				output.add("@NodeStyle\tTaxonomy\tCircle");
				output.add("@ColorTable\tFews8\tWhite-Green");
				output.add("@ColorEdits");	
				for(int key:table.getClass2Counts("Taxonomy").keySet()){
					output.add("TAX\t"+key+"\t"+table.getClass2Counts("Taxonomy").get(key)[0]);
				}
				output.add("END_OF_DATA_TABLE");
				output.add("#SampleID\t" + "@Source\t"+"Size");
				output.add(fileName+"\t"+file.getFileName()+"\t"+table.getSampleSizes()[0]);
				Path path = Paths.get(outDir+"/MeganSummaries/"+fileName+".megan");
				Files.write(path, output, Charset.forName("UTF-8"));
			}else{
				warning.log(Level.SEVERE,"File does not exist: "+fileName);
			}
		} catch (IOException e) {
			warning.log(Level.SEVERE,"Cannot Write or Read file",e);
		}
	}

}
