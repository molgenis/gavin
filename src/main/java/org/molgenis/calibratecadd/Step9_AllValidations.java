package org.molgenis.calibratecadd;

import java.io.File;
import java.util.Arrays;
import java.util.List;


public class Step9_AllValidations
{
	public static void main(String[] args) throws Exception
	{
		String outFile = "/Users/jvelde/toolcomparison_output.r";
		String geneSumm = "/Users/jvelde/github/maven/molgenis-data-cadd/data/predictions/GAVIN_calibrations_r0.1.tsv";
		List<String> datasets = Arrays.asList(new String[]{"ClinVarNew", "MutationTaster2", "UMCG_Onco", "UMCG_Various", "VariBenchTest", "VariBenchTraining"});
		List<String> tools = Arrays.asList(new String[]{"GAVIN", "PONP2", "CADD", "PROVEAN", "SIFT", "PolyPhen2", "MSC", "Condel"});

		if(new File(outFile).exists())
		{
			throw new Exception("output file already exists: " + outFile);
		}
		new File(outFile).createNewFile();

		for(String dataset: datasets)
		{
			for(String tool: tools)
			{
				new Step9_Validation(geneSumm, "/Users/jvelde/github/maven/molgenis-data-cadd/data/goldstandards/" + dataset, tool, outFile);
			}
		}

	}
	
}
