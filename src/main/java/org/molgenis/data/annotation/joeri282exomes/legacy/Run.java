package org.molgenis.data.annotation.joeri282exomes.legacy;

import java.io.File;

public class Run
{

	/**
	 * Check all arguments and start the tool
	 * 
	 * @param args
	 * @throws Exception
	 */
	public static void main(String[] args) throws Exception
	{
		
		if (!(args.length == 2))
		{
			throw new Exception(
					"Must supply 2 arguments");
		}
		
		File vcfFile = new File(args[0]);
		if (!vcfFile.isFile())
		{
			throw new Exception("Input VCF file does not exist or directory: " + vcfFile.getAbsolutePath());
		}
		
		File patientGroups = new File(args[1]);
		if (!patientGroups.isFile())
		{
			throw new Exception("Patient groups file does not exist or directory: " + patientGroups.getAbsolutePath());
		}

		DifferentialExomes de = new DifferentialExomes(vcfFile, patientGroups);
		de.start();

	}

}
