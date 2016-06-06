package org.molgenis.calibratecadd.support;

import org.molgenis.calibratecadd.Step9_Validation;

import java.io.File;
import java.util.List;


public class BootStrappingAnalysisRunner
{
	public static void main(String[] args) throws Exception
	{
		String fullSet = "/Users/joeri/Desktop/old/Desktop/clinvarcadd/combined_datasets_for_external_scoring/cat_all_vcfs.txt";
		String gavin = "/Users/joeri/github/gavin/data/predictions/GAVIN_calibrations_r0.1.tsv";
		String outFile = "/Users/joeri/github/gavin/data/other/performancebootstrap_output_usedinpaper.r";

		int iterations = 10000;

		if(new File(outFile).exists())
		{
			throw new Exception("output file already exists: " + outFile);
		}
		new File(outFile).createNewFile();

		BootStrappingAnalysis ba = new BootStrappingAnalysis(fullSet, gavin, outFile, Step9_Validation.ToolNames.GAVIN);

		for(int i = 0 ; i < iterations; i ++)
		{
			ba.getStatsOnSet(ba.randomSubset(100, 100, BootStrappingVariant.Label.C1_C2));
			ba.getStatsOnSet(ba.randomSubset(100, 100, BootStrappingVariant.Label.C3));
			ba.getStatsOnSet(ba.randomSubset(100, 100, BootStrappingVariant.Label.C4));
		}

		ba = new BootStrappingAnalysis(fullSet, gavin, outFile, Step9_Validation.ToolNames.GAVINnocal);

		for(int i = 0 ; i < iterations; i ++)
		{
			ba.getStatsOnSet(ba.randomSubset(100, 100, BootStrappingVariant.Label.C1_C2));
			ba.getStatsOnSet(ba.randomSubset(100, 100, BootStrappingVariant.Label.C3));
			ba.getStatsOnSet(ba.randomSubset(100, 100, BootStrappingVariant.Label.C4));
		}
	}
}
