package org.molgenis.calibratecadd;

import java.io.File;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.util.Arrays;
import java.util.List;


public class Step9_AllValidations
{
	public static void main(String[] args) throws Exception
	{
		//paper submitted: r0.1 (using clinvar: 5 nov 2015). after revision: r0.2, using the same clinvar data but snpeff-annotated ExAC
		//latest: r0.3 (using clinvar: 1 sep 2016)
		String version = "r0.2";
		String outFile = "/Users/joeri/github/gavin/data/other/step9_panels_out_"+version+".R";
		String path = "/Users/joeri/github/gavin/data/predictions";
		List<String> datasets = Arrays.asList(new String[]{"Renal", "Pulmonary", "Ophthalmologic", "Oncologic", "Obstetric", "NotInCGD", "Neurologic", "Musculoskeletal", "Hematologic", "Genitourinary", "Gastrointestinal", "Endocrine", "Dermatologic", "Dental", "Craniofacial", "Cardiovascular", "Biochemical", "Audiologic_Otolaryngologic", "Allergy_Immunology_Infectious"});
		List<String> tools = Arrays.asList(new String[]{"GAVIN", "CADD", "MSC_ClinVar95CI", "MSC_HGMD99CI", "PROVEAN", "SIFT", "PolyPhen2", "Condel", "PONP2", "PredictSNP2", "FATHMM", "GWAVA", "FunSeq", "DANN"});

		if(new File(outFile).exists())
		{
			throw new Exception("output file already exists: " + outFile);
		}
		new File(outFile).createNewFile();
		Files.write(Paths.get(outFile), "df <- data.frame()\n".getBytes(), StandardOpenOption.APPEND);

		for(String dataset: datasets)
		{
			for(String tool: tools)
			{
				new Step9_Validation(path, "/Users/joeri/github/gavin/data/goldstandards/cgdpanels/" + dataset, tool, outFile, version);
			}
		}

	}
	
}
