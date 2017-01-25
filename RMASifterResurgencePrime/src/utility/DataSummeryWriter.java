package utility;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;

import megan.chart.ChartColorManager;
import megan.core.DataTable;
import megan.core.Director;
import megan.core.Document;
import megan.core.MeganFile;
import megan.core.SampleAttributeTable;
import megan.core.SyncArchiveAndDataTable;
import megan.data.IConnector;
import megan.rma6.RMA6File;

public class DataSummeryWriter {

	public static void main(String[] args) throws IOException {
		// TODO Auto-generated method stub
		Document doc = new Document();
		SampleAttributeTable sampleAttributeTable = new SampleAttributeTable();
		DataTable table = new DataTable();
		MeganFile file = new MeganFile();
		file.setFileFromExistingFile("/Users/huebler/Desktop/Quid_Data/RMA6/quid.AR.collapsed.rma6", true);
		file.setReadOnly(true);
		
//		System.out.println(doc.getNumberOfSamples());
//		System.out.println(file.getDataConnector(true).getNumberOfMatches());
		IConnector connector = file.getDataConnector();
        SyncArchiveAndDataTable.syncArchive2Summary(file.getFileName(), connector, table, sampleAttributeTable);
        doc.parseParameterString(table.getParameters());
//        System.out.println(doc.getParameterString());
//		System.out.println(table.getColorTable());
		ChartColorManager chartColorManager = new  ChartColorManager(null);
       System.out.println("@Creator\t"+ table.getCreator());
       System.out.println("@CreationDate\t"+ table.getCreationDate());
       System.out.println("@ContentType\t"+table.getContentType());
       System.out.println("@Names\t"+table.getSampleNames()[0]);
       System.out.println("@BlastMode\t"+table.getBlastMode());
       System.out.println("@Uids\t"+table.getSampleUIds()[0]);
       System.out.println("@Sizes\t"+table.getSampleSizes()[0]);
       System.out.println("@TotalReads\t"+table.getTotalReads());
       System.out.println("@AdditionalReads\t"+table.getAdditionalReads());
       System.out.println("@Collapse\t"+"Taxonomy\t"+"-1\t");//table.getCollapsed("Taxonomy"));//what are you
       System.out.println("@Algorithm\t"+table.getAlgorithm("Taxonomy"));
       System.out.println("@NodeStyle\tTaxonomy\tCircle");
       System.out.println("@ColorTable\tFews8\tWhite-Green");
       System.out.print("@ColorEdits");	
       for(int key:table.getClass2Counts("Taxonomy").keySet()){
    	   System.out.println("TAX\t"+key+"\t"+table.getClass2Counts("Taxonomy").get(key)[0]);
       }
       System.out.println("END_OF_DATA_TABLE");
       System.out.println("#SampleID\t" + "@Source\t"+"Size");
       System.out.println(file.getName()+"\t"+file.getFileName()+"\t"+table.getSampleSizes()[0]);

	}

}
