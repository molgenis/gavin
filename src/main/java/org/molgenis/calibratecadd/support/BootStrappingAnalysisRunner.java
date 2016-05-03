package org.molgenis.calibratecadd.support;

import java.util.List;


public class BootStrappingAnalysisRunner
{
	public static void main(String[] args) throws Exception
	{
		String fullSet = "/Users/jvelde/Desktop/clinvarcadd/combined_datasets_for_external_scoring/cat_all_vcfs.txt";
		String gavin = "/Users/jvelde/github/maven/molgenis-data-cadd/data/predictions/GAVIN_calibrations_r0.1.tsv";
		String outFile = "/Users/jvelde/bootstrapresults.r";
		
		int iterations = 10000;
		int sampleSize = 10000;
		
		BootStrappingAnalysis ba = new BootStrappingAnalysis(fullSet, gavin, outFile);
//		ba.getStatsOnFullSet();
		for(int i = 0 ; i < iterations; i ++)
		{
			List<BootStrappingVariant> set = ba.randomSubset(sampleSize);
			List<BootStrappingVariant> dsSet = ba.downSampleToUniformCalibPercDistr(set);
			ba.getStatsOnSet(dsSet);
		}
			
		
	}
}
