package org.molgenis.caddtlmapping.binom;

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
		if (!(args.length == 4 || args.length == 5))
		{
			throw new Exception(
					"Must supply at least 5 arguments: patient VCF (tabix indexed), ExAC file, CADD folder, inheritance mode, optional: sample identifier list\n"
							+ "Typical input prep:\n"
							+ "java -Xmx2g -jar snpEff.jar hg19 -v -canon -ud 0 1504.vcf > 1504_snpeff.vcf\n"
							+ "vcf-sort 1504_snpeff.vcf | bgzip > 1504_snpeff_sorted.vcf.gz\n"
							+ "tabix -p vcf 1504_snpeff_sorted.vcf.gz\n" + "java -Xmx2g -jar CADDTLMapping .. etc\n");
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
		if (!caddFile.isDirectory())
		{
			throw new Exception("Input CADD folder does not exist or directory: " + caddFile.getAbsolutePath());
		}

		String inheritance = args[3];
		if (!(inheritance.equals("additive") || inheritance.equals("dominant") || inheritance.equals("recessive")))
		{
			throw new Exception("Input inheritance model must be 'additive', 'dominant' or 'recessive', not "
					+ inheritance);
		}

		File patientSampleIdsFile = null;
		if (args.length == 5)
		{
			patientSampleIdsFile = new File(args[4]);
			if (!patientSampleIdsFile.isFile())
			{
				throw new Exception("Input patient sample IDs file does not exist or directory: "
						+ patientSampleIdsFile.getAbsolutePath());
			}
		}

		System.out.println("Arguments OK !\nStarting..");

		CADDTLMapping cm = new CADDTLMapping(vcfFile, exacFile, caddFile, inheritance, patientSampleIdsFile);
		cm.start();

		System.out.println("All done!");
	}

}
