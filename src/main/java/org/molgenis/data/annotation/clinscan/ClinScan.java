package org.molgenis.data.annotation.clinscan;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Map;
import java.util.Scanner;
import java.util.Set;
import java.util.stream.Stream;

import org.molgenis.data.Entity;
import org.molgenis.data.annotation.RepositoryAnnotator;
import org.molgenis.data.annotation.cmd.CommandLineAnnotatorConfig;
import org.molgenis.data.vcf.VcfRepository;
import org.molgenis.util.ApplicationContextProvider;
import org.springframework.context.ApplicationContext;
import org.springframework.context.annotation.AnnotationConfigApplicationContext;
import org.springframework.stereotype.Component;

@Component
public class ClinScan
{

	private String gtcMessage = null;

	public void go(File vcfFile, File exacFile) throws Exception
	{
		VcfRepository vcf = new VcfRepository(vcfFile, "vcf");
		ApplicationContext applicationContext = ApplicationContextProvider.getApplicationContext();
		Map<String, RepositoryAnnotator> annotators = applicationContext.getBeansOfType(RepositoryAnnotator.class);
		RepositoryAnnotator exacAnnotator = annotators.get("exac");

		exacAnnotator.getCmdLineAnnotatorSettingsConfigurer().addSettings(exacFile.getAbsolutePath());
		// annotate(exacAnnotator, inputVcfFile, outputVCFFile, attrNames);

		Iterator<Entity> it = exacAnnotator.annotate(vcf);

		while (it.hasNext())
		{
			Entity record = it.next();

			String filter = record.getString("FILTER");

			if (!filter.equals("PASS"))
			{
				continue;
			}

			String chr = record.getString("#CHROM");
			String pos = record.getString("POS");
			String ref = record.getString("REF");
			String altStr = record.getString("ALT");
			String exac_af_STR = record.get("EXAC_AF") == null ? null : record.get("EXAC_AF").toString();
			String exac_ac_hom_STR = record.get("EXAC_AC_HOM") == null ? null : record.get("EXAC_AC_HOM").toString();
			String exac_ac_het_STR = record.get("EXAC_AC_HET") == null ? null : record.get("EXAC_AC_HET").toString();
			String[] multiAnn = record.getString("ANN").split(",");

			String[] altsplit = altStr.split(",", -1);

			String[] exac_af_split = new String[altsplit.length];
			if (exac_af_STR != null)
			{
				exac_af_split = exac_af_STR.split(",", -1);
			}

			String[] exac_ac_hom_split = new String[altsplit.length];
			if (exac_ac_hom_STR != null)
			{
				exac_ac_hom_split = exac_ac_hom_STR.split(",", -1);
			}

			String[] exac_ac_het_split = new String[altsplit.length];
			if (exac_ac_het_STR != null)
			{
				exac_ac_het_split = exac_ac_het_STR.split(",", -1);
			}

			/**
			 * Iterate over alternatives, if applicable multi allelic example: 1:1148100-1148100
			 */
			for (int i = 0; i < altsplit.length; i++)
			{
				String alt = altsplit[i];

				if (exac_af_STR != null && !exac_af_split[i].equals("."))
				{
					Double exac_af = Double.parseDouble(exac_af_split[i]);
					if (exac_af > 0.05)
					{
						continue;
					}
				}

				String ann = multiAnn[i];
				String[] annSplit = ann.split("\\|", -1);
				String impact = annSplit[2];
				if (impact.equals("MODIFIER") || impact.equals("LOW"))
				{
					continue;
				}

				String effect = annSplit[1];

				String cDNA = annSplit[9];
				String aaChange = annSplit[10];

				String gene = annSplit[3];
				if (gene.isEmpty())
				{
					throw new Exception("reminder: gene can be empty??");
				}

				// TODO: add CADD scores

				int[] patient_GTC = countGTC(record, i);
				if (patient_GTC == null)
				{
					continue;
				}

				String ExAC_AF = exac_af_split[i] == null ? "0" : exac_af_split[i];
				int ExAC_AC_HOM = exac_ac_hom_split[i] == null ? 0 : Integer.parseInt(exac_ac_hom_split[i]);
				int ExAC_AC_HET = exac_ac_het_split[i] == null ? 0 : Integer.parseInt(exac_ac_het_split[i]);

				// if not actually seen in patients, skip it...
				if (patient_GTC[1] == 0 && patient_GTC[2] == 0)
				{
					continue;
				}

				// if the number of people in exac with het or homalt exceed our patients with this variant.. it's hard
				// to believe it :-\ exac has late-onset, common diseases, see: http://exac.broadinstitute.org/about
				// if (ExAC_AC_HET > patient_GTC[1] || ExAC_AC_HOM > patient_GTC[2])
				if (ExAC_AC_HET > 10 || ExAC_AC_HOM > 10)
				{
					continue;
				}

				String variantInfo = chr + ":" + pos + "-" + pos + ", " + ref + "/" + alt + ", " + cDNA + ", "
						+ aaChange + ", " + gene + ", effect: " + effect + ", impact: " + impact
						+ ", ExAC [allelefreq=" + ExAC_AF + ", hets=" + ExAC_AC_HET + ", homalts=" + ExAC_AC_HOM
						+ "], patients [homrefs=" + patient_GTC[0] + ", hets=" + patient_GTC[1] + ", homalts="
						+ patient_GTC[2] + "], details: [" + this.gtcMessage + "]";

				System.out.println("candidate: " + variantInfo);

			}

		}
	}

	public int[] countGTC(Entity record, int altIndex)
	{
		this.gtcMessage = "";

		// because alt index = 0 for the first alt, we add 1
		altIndex = altIndex + 1;

		// for a particular ref-alt combination:
		// [homref, het, homalt]
		// can only do for ref/alt-index combinations, so e.g. 0/0, 0|2 or 2/2. print warning on 1/2, 3|2, etc.
		// warn if this happens
		// also count for "other people" in [3][4] and [5] as control reference
		int[] gtc = new int[]
		{ 0, 0, 0 };

		Iterable<Entity> sampleEntities = record.getEntities(VcfRepository.SAMPLES);
		for (Entity sample : sampleEntities)
		{
			String sampleName = sample.get("ORIGINAL_NAME").toString();

			String genotype = sample.get("GT").toString();

			if (genotype.equals("./."))
			{
				continue;
			}

			// quality filter: we want depth X or more
			int depthOfCoverage = Integer.parseInt(sample.get("DP").toString());
			if (depthOfCoverage < 10)
			{
				continue;
			}

			if (genotype.equals("0/0") || genotype.equals("0|0"))
			{

				gtc[0]++;
			}
			else if (genotype.equals("0/" + altIndex) || genotype.equals(altIndex + "/0")
					|| genotype.equals("0|" + altIndex) || genotype.equals(altIndex + "|0"))
			{
				gtc[1]++;
				this.gtcMessage += "pathet:" + sampleName + ",dp:" + depthOfCoverage + ",gt:" + genotype + " ";
			}
			else if (genotype.equals(altIndex + "/" + altIndex) || genotype.equals(altIndex + "|" + altIndex))
			{
				gtc[2]++;
				this.gtcMessage += "pathom:" + sampleName + ",dp:" + depthOfCoverage + ",gt:" + genotype + " ";
			}
			else if (genotype.contains(altIndex + ""))
			{
				System.out.println("WARNING: genotype " + genotype + " not counted for altindex " + altIndex + " for "
						+ record.getString("#CHROM") + ":" + record.getString("POS")
						+ " because it's not a ref-alt combination!");
			}

		}

		return gtc;

	}

	public static void main(String[] args) throws Exception
	{
		// configureLogging();

		// See http://stackoverflow.com/questions/4787719/spring-console-application-configured-using-annotations
		AnnotationConfigApplicationContext ctx = new AnnotationConfigApplicationContext("org.molgenis.data.annotation");
		ctx.register(CommandLineAnnotatorConfig.class);
		ClinScan main = ctx.getBean(ClinScan.class);
		main.run(args);
		ctx.close();
	}

	public void run(String[] args) throws Exception
	{
		if (!(args.length == 2))
		{
			throw new Exception("Must supply 2 arguments");
		}

		File vcfFile = new File(args[0]);
		if (!vcfFile.isFile())
		{
			throw new Exception("Input VCF file does not exist or directory: " + vcfFile.getAbsolutePath());
		}

		File exacFile = new File(args[1]);
		if (!exacFile.isFile())
		{
			throw new Exception("Exac file does not exists at " + exacFile);
		}

		ClinScan cf = new ClinScan();
		cf.go(vcfFile, exacFile);

	}

}
