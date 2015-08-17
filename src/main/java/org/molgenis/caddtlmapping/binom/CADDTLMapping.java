package org.molgenis.caddtlmapping.binom;

import java.io.File;

/**
 * Fresh implementation based on binomial test of over represented genotypes for variants within a certain range of CADD
 * scores. Does not need MAF cutoff to find candidates because we want to use the underlying genotype data to do the
 * statistical test. Takes a folder of CADD files, so we can download multiple and check them all because we want to
 * assess as many variants as possible.
 * 
 * @author jvelde
 *
 */
public class CADDTLMapping
{

	public CADDTLMapping(File vcfFile, File exacFile, File caddFolder, String inheritance, File patientSampleIdsFile)
	{

	}

	public void start() throws Exception
	{

	}
}
