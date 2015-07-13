package org.molgenis.caddtlmapping;

import java.io.File;

public class Run
{

	/**
	 * Check all arguments and start the tool
	 * @param args
	 * @throws Exception
	 */
	public static void main(String[] args) throws Exception
	{
		if (!(args.length == 6 || args.length == 7))
		{
			throw new Exception(
					"Must supply at least 6 arguments: patient VCF, ExAC, CADD, output, MAF, inheritance, sampleIds");
		}

		File vcfFile = new File(args[0]);
		if (!vcfFile.isFile())
		{
			throw new Exception("Input VCF file does not exist or directory: " + vcfFile.getAbsolutePath());
		}

		File exacFile = new File(args[1]);
		if (!exacFile.isFile())
		{
			throw new Exception("Input ExAC file does not exist or directory: " + exacFile.getAbsolutePath());
		}

		File caddFile = new File(args[2]);
		if (!caddFile.isFile())
		{
			throw new Exception("Input CADD file does not exist or directory: " + caddFile.getAbsolutePath());
		}

		File outputFile = new File(args[3]);
		if (outputFile.isFile())
		{
			System.out.println("WARNING: output file already exists: " + outputFile.getAbsolutePath());
		}

		String MAF_input = args[4];
		double MAF = -1;
		try
		{
			MAF = Double.parseDouble(MAF_input);
		}
		catch (NumberFormatException e)
		{
			throw new Exception("Input MAF not a valid double value: " + MAF_input);
		}
		if (MAF < 0 || MAF > 1)
		{
			throw new Exception("Input MAF must be between 0 and 1:" + MAF_input);
		}

		String inheritance = args[5];
		if (!(inheritance.equals("additive") || inheritance.equals("dominant") || inheritance.equals("recessive")))
		{
			throw new Exception("Input inheritance model must be 'additive', 'dominant' or 'recessive', not "
					+ inheritance);
		}

		File patientSampleIdsFile = null;
		if (args.length == 7)
		{
			patientSampleIdsFile = new File(args[6]);
			if (!patientSampleIdsFile.isFile())
			{
				throw new Exception("Input patient sample IDs file does not exist or directory: "
						+ patientSampleIdsFile.getAbsolutePath());
			}
		}

		System.out.println("Arguments OK !\nStarting..");

		CADDTLMapping cm = new CADDTLMapping(vcfFile, exacFile, caddFile, outputFile, MAF, inheritance,
				patientSampleIdsFile);
		cm.start();
		
		System.out.println("All done!");
	}

}
