package org.molgenis.calibratecadd;

import java.io.File;
import java.util.Arrays;
import java.util.List;


public class Step9_AllValidations
{
	public static void main(String[] args) throws Exception
	{
		String outFile = "/Users/joeri/step9_all_out.R";
		String path = "/Users/joeri/github/gavin/data/predictions";
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
				new Step9_Validation(path, "/Users/joeri/github/gavin/data/goldstandards/" + dataset, tool, outFile);
			}
		}

	}
	
}
