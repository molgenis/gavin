package org.molgenis.data.annotation.joeri282exomes;

import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Scanner;
import java.util.Set;

import org.molgenis.data.Entity;
import org.molgenis.data.annotation.cmd.CommandLineAnnotatorConfig;
import org.molgenis.data.vcf.VcfRepository;
import org.springframework.context.annotation.AnnotationConfigApplicationContext;
import org.springframework.stereotype.Component;

@Component
public class CalibratedExomeAnalysis
{

	HashMap<String, ArrayList<String>> patientgroupToGenes = new HashMap<String, ArrayList<String>>();

	private HashMap<String, String> sampleToGroup;

	private String gtcMessage = null;

	public void go(File vcfFile, File patientGroups) throws Exception
	{

		Set<String> groups = new HashSet<String>();

		HashMap<String, String> sampleToGroup = new HashMap<String, String>();
		Scanner s = new Scanner(patientGroups);
		String line = null;
		while (s.hasNextLine())
		{
			line = s.nextLine();
			String[] lineSplit = line.split("\t", -1);
			sampleToGroup.put(lineSplit[0], lineSplit[1]);
			groups.add(lineSplit[1]);
		}
		s.close();

		this.sampleToGroup = sampleToGroup;

		VcfRepository vcf = new VcfRepository(vcfFile, "vcf");
		Iterator<Entity> it = vcf.iterator();

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

			//TODO: multiple genes!?!? check MTOR / ANGPTL7
			//e.g.
			//1       11249803        rs141187277     T       C       1634.44 PASS    AC=2;AF=0.003546;AN=564;BaseQRankSum=4.12;ClippingRankSum=0.209;DB;DP=10525;FS=0.0;GQ_MEAN=98.68;GQ_STDDEV=45.19;InbreedingCoeff=-0.0036;MLEAC=2;MLEAF=0.003546;MQ=60.0;MQ0=0;MQRankSum=1.01;NCC=0;QD=17.96;ReadPosRankSum=0.377;SOR=0.622;VQSLOD=18.2;culprit=MQ;ANN=C|missense_variant|MODERATE|ANGPTL7|ANGPTL7|transcript|NM_021146.3|Coding|1/5|c.167T>C|p.Val56Ala|458/2290|167/1041|56/346||,C|intron_variant|MODIFIER|MTOR|MTOR|transcript|NM_004958.3|Coding|28/57|c.4253+9512A>G||||||;GoNL_GTC=.,.,.;GoNL_AF=.;EXAC_AF=1.318E-4;EXAC_AC_HOM=0;EXAC_AC_HET=16;CADD=3.259127;CADD_SCALED=22.8;
			
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


			}

		}
	}

	public static void main(String[] args) throws Exception
	{
		// configureLogging();

		// See http://stackoverflow.com/questions/4787719/spring-console-application-configured-using-annotations
		AnnotationConfigApplicationContext ctx = new AnnotationConfigApplicationContext("org.molgenis.data.annotation");
		ctx.register(CommandLineAnnotatorConfig.class);
		CalibratedExomeAnalysis main = ctx.getBean(CalibratedExomeAnalysis.class);
		main.run(args);
		ctx.close();
	}

	public void run(String[] args) throws Exception
	{
		//TODO: first pass, print list of variants that miss a CADD score
		//we can do SNVs easy but not indels and other 'non predictable' changes
		//send off to CADD webservice, and give the output so we use this in second pass
		
		if (!(args.length == 2))
		{
			throw new Exception("Must supply 2 arguments");
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

		CalibratedExomeAnalysis cf = new CalibratedExomeAnalysis();
		cf.go(vcfFile, patientGroups);

	}

}
