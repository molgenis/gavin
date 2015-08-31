package org.molgenis.data.annotation.joeri282exomes;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Map;
import java.util.Scanner;
import java.util.Set;

import org.molgenis.data.Entity;
import org.molgenis.data.annotation.RepositoryAnnotator;
import org.molgenis.data.annotation.cmd.CommandLineAnnotatorConfig;
import org.molgenis.data.vcf.VcfRepository;
import org.molgenis.util.ApplicationContextProvider;
import org.springframework.context.ApplicationContext;
import org.springframework.context.annotation.AnnotationConfigApplicationContext;
import org.springframework.stereotype.Component;

@Component
public class ClinicalFilters
{
	/**
	 * Candidates
	 */
	// source: CGD, matching 'stickler syndrome'
	ArrayList<String> sticklerGenes = new ArrayList<String>(Arrays.asList(new String[]
	{ "COL11A1", "COL11A2", "COL9A1", "COL9A2", "COL2A1" }));

	// source: CGD, matching 'glycogen storage disease'
	ArrayList<String> gsdGenes = new ArrayList<String>(Arrays.asList(new String[]
	{ "AGL", "ALDOA", "ENO3", "G6PC", "GAA", "GBE1", "GYG1", "GYS1", "GYS2", "LDHA", "PFKM", "PGAM2", "PHKA1", "PHKA2",
			"PHKB", "PHKG2", "PRKAG2", "PYGL", "PYGM", "SLC2A2", "SLC37A4" }));

	// source: http://www.ncbi.nlm.nih.gov/books/NBK121988/ - 6 of which in CGD as well
	ArrayList<String> ironaccGenes = new ArrayList<String>(Arrays.asList(new String[]
	{ "PANK2", "PLA2G6", "C19orf12", "FA2H", "ATP13A2", "WDR45", "COASY", "FTL", "CP", "DCAF17" }));

	// source: CGD, matching 'spinal muscular atrophy'
	ArrayList<String> smaGenes = new ArrayList<String>(
			Arrays.asList(new String[]
			{ "ASAH1", "ATP7A", "CHCHD10", "DNAJB2", "DYNC1H1", "IGHMBP2", "PLEKHG5", "SMN1", "SMN2", "TRPV4", "UBA1",
					"VAPB" }));

	private HashMap<String, String> sampleToGroup;
	private Set<String> groups;

	public void go(File vcfFile, File patientGroups, File exacFile) throws Exception
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
		this.groups = groups;

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
			String[] multiAnn = record.getString("INFO_ANN").split(",");

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

				String gene = annSplit[3];
				if (gene.isEmpty())
				{
					throw new Exception("reminder: gene can be empty??");
				}

				// TODO: print genotype counts
				// TODO: print CADD scores

				int[] patient_GTC = countGTC(record, i);

				String ExAC_AC_HOM = exac_ac_hom_split[i] == null ? "0" : exac_ac_hom_split[i];
				String ExAC_AF = exac_af_split[i] == null ? "0" : exac_af_split[i];

				String variantInfo = chr + ":" + pos + "-" + pos + ", " + ref + "/" + alt + ", " + gene + ", ExAC_AF="
						+ ExAC_AF + " (ExAC hom=" + ExAC_AC_HOM + ", Patient hom=" + patient_GTC[2] + "), impact: "
						+ impact + "[pat homref=" + patient_GTC[0] + " het=" + patient_GTC[1] + "]";

				if (sticklerGenes.contains(gene))
				{
					System.out.println("stickler candidate? " + variantInfo);
				}

				if (gsdGenes.contains(gene))
				{
					System.out.println("gsd candidate? " + variantInfo);
				}

				if (ironaccGenes.contains(gene))
				{
					System.out.println("ironacc candidate? " + variantInfo);
				}

				if (smaGenes.contains(gene))
				{
					System.out.println("sma candidate? " + variantInfo);
				}

			}

		}
	}

	public int[] countGTC(Entity record, int altIndex)
	{
		// because alt index = 0 for the first alt, we add 1
		altIndex = altIndex + 1;

		// for a particular ref-alt combination:
		// [homref, het, homalt]
		// can only do for ref/alt-index combinations, so e.g. 0/0, 0|2 or 2/2. print warning on 1/2, 3|2, etc.
		// error if this happens
		int[] gtc = new int[]
		{ 0, 0, 0 };

		Iterable<Entity> sampleEntities = record.getEntities(VcfRepository.SAMPLES);
		for (Entity sample : sampleEntities)
		{

			String genotype = sample.get("GT").toString();

			if (genotype.equals("./."))
			{
				continue;
			}

			// TODO: quality filter
			
			// TODO: filter on patient_groups for GTC (scan for each?)

			if (genotype.equals("0/0") || genotype.equals("0|0"))
			{
				gtc[0]++;
			}
			else if (genotype.equals("0/" + altIndex) || genotype.equals("0|" + altIndex)
					|| genotype.equals(altIndex + "|0"))
			{
				gtc[1]++;
			}
			else if (genotype.equals(altIndex + "/" + altIndex) || genotype.equals(altIndex + "|" + altIndex))
			{
				gtc[2]++;
			}
			else if (genotype.contains(altIndex + ""))
			{
				System.out.println("WARNING: genotype " + genotype + " not counted for altindex " + altIndex + " for "
						+ record.getString("#CHROM") + ":" + record.getString("POS") + " because it's not a ref-alt combination!");
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
		ClinicalFilters main = ctx.getBean(ClinicalFilters.class);
		main.run(args);
		ctx.close();
	}

	public void run(String[] args) throws Exception
	{
		if (!(args.length == 3))
		{
			throw new Exception("Must supply 3 arguments");
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

		File exacFile = new File(args[2]);
		if (!exacFile.isFile())
		{
			throw new Exception("Exac file does not exists at " + exacFile);
		}

		ClinicalFilters cf = new ClinicalFilters();
		cf.go(vcfFile, patientGroups, exacFile);

	}

}
