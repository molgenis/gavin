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
		String outFile = "/Users/joeri/github/gavin/data/other/step9_panels_out.R";
		String path = "/Users/joeri/github/gavin/data/predictions";
	//	List<String> datasets = Arrays.asList(new String[]{"ClinVarNew", "MutationTaster2", "UMCG_Onco", "UMCG_Various", "VariBenchTest", "VariBenchTraining"});
		List<String> datasets = Arrays.asList(new String[]{"Renal", "Pulmonary", "Ophthalmologic", "Oncologic", "Obstetric", "NotInCGD", "Neurologic", "Musculoskeletal", "Hematologic", "Genitourinary", "Gastrointestinal", "Endocrine", "Dermatologic", "Dental", "Craniofacial", "Cardiovascular", "Biochemical", "Audiologic_Otolaryngologic", "Allergy_Immunology_Infectious"});
		List<String> tools = Arrays.asList(new String[]{"GAVIN", "CADD", "MSC_ClinVar95CI", "MSC_HGMD99CI", "PROVEAN", "SIFT", "PolyPhen2", "Condel", "PONP2"});

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
//				new Step9_Validation(path, "/Users/joeri/github/gavin/data/goldstandards/" + dataset, tool, outFile);
				new Step9_Validation(path, "/Users/joeri/github/gavin/data/goldstandards/cgdpanels/" + dataset, tool, outFile);
			}
		}

	}
	
}
