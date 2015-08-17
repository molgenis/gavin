package org.molgenis.caddtlmapping.mannwhitney;

import java.io.File;
import java.io.PrintWriter;
import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Date;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;

import org.apache.commons.math3.stat.inference.MannWhitneyUTest;

// "CADD-TL mapping"
// test 3 inheritance models: recessive, additive, dominant
// test 3 frequency models: singleton (0%), rare (<1%), common (<5%)
//
// requires patient VCF to be indexed TODO: bgzip myfile.vcf; tabix -p vcf myfile.vcf.gz; 
//
// simple mechanics:
// 1. get set of candidate variants (+genotypes) from patients by filtering on MAF, impact, etc
// 2. from within the same sequence feature (ie. gene), get comparable variants from the population
// meaning: same filter criteria. get many, but not more than necessary while staying close to patient variants
// ---> QUESTION: okay if the same, with different genotypes???
// 3. get CADD scores and multiply by inheritance model
// (homalt always x2, but rec: 0x het, add: 1x het, dom: 2x het)
// 4. test and plot!
//
// CURRENT LIMITATIONS/TODO:
// - Exome only, later on use 1000G+GoNL data for whole-genome scan (need to include feature locations for enhancers, TFBS, lncRNA etc ??)
// - Ignoring multi-allelic variants for the moment as they are tricky to parse
// - For developing, only use subset
// - handle hemizygous on X / Y ??
// - see TODO inline
// - see ftp://ftp.broadinstitute.org/pub/ExAC_release/release0.3/README.known_issues
// - cleanup dependencies in POM file
//
public class CADDTLMapping
{
	private File vcfFile;
	private File exacFile;
	private File caddFile; // TODO directory of diffent CADD files to maximize succes
	private File outputFile;
	private double mafThreshold;
	private String inheritance;
	private File patientSampleIdsFile;

	public CADDTLMapping(File vcfFile, File exacFile, File caddFile, double mafThreshold, String inheritance,
			File patientSampleIdsFile)
	{
		super();
		this.vcfFile = vcfFile;
		this.exacFile = exacFile;
		this.caddFile = caddFile;
		DateFormat dateFormat = new SimpleDateFormat("yyyy-MM-dd-HH-mm-ss");
		Date date = new Date();
		this.outputFile = new File("caddtlmappingresults_" + dateFormat.format(date) + "_" + inheritance + "_"
				+ mafThreshold + ".tsv");
		this.mafThreshold = mafThreshold;
		this.inheritance = inheritance;
		this.patientSampleIdsFile = patientSampleIdsFile;
	}

	/**
	 * Central function of the tool
	 * 
	 * @throws Exception
	 */
	public void start() throws Exception
	{

		Helper h = new Helper(vcfFile, exacFile, caddFile);

		// list of sample identifers within the VCF file that represent the patients
		// TODO: check VCF file and give warnings when there are identifiers not mapping etc
		ArrayList<String> patientSampleIdList = null;
		if (patientSampleIdsFile != null)
		{
			System.out.println("Loading patient sampleID list..");
			patientSampleIdList = h.loadPatientSampleIdList(patientSampleIdsFile);
		}

		System.out.println("First pass to tag potentially pathogenic variants and store feature locations..");
		LinkedHashMap<String, String> sequenceFeatureLocations = new LinkedHashMap<String, String>();
		LinkedHashMap<String, List<String>> candidates = h.getCandidateVariantsPerGene(vcfFile, mafThreshold,
				sequenceFeatureLocations);
		System.out.println("A total of " + candidates.size() + " sequence features have 1 or more tagged variants.");

		System.out.println("Second pass to perform CADDTL mapping on each sequence feature..");
		HashMap<String, Double> lodScores = new HashMap<String, Double>();
		HashMap<String, Double> pvals = new HashMap<String, Double>();

		for (String sequenceFeature : candidates.keySet())
		{
			System.out.println("Now investigating " + sequenceFeature + " (" + candidates.get(sequenceFeature).size()
					+ " tagged variants)");

			double[] patientCaddScores = h.getPatientCaddScores(candidates.get(sequenceFeature), inheritance,
					patientSampleIdList);
			double[] populationCaddScores = h.getPopulationCaddScores(candidates.get(sequenceFeature), inheritance);

			double lod = 0;
			double pval = 1;
			if (patientCaddScores.length > 0 && populationCaddScores.length > 0)
			{
				MannWhitneyUTest mwt = new MannWhitneyUTest();
				pval = mwt.mannWhitneyUTest(patientCaddScores, populationCaddScores);
				pval = pval != 0 ? pval : 0.0000000001;
				lod = -Math.log10(pval);
			}

			lodScores.put(sequenceFeature, lod);
			pvals.put(sequenceFeature, pval);

			System.out.println("Evaluated " + patientCaddScores.length + " patient scores vs. "
					+ populationCaddScores.length + " population scores resulting in a LOD score of " + lod);
		}

		System.out.println("Processing the results..");
		LinkedHashMap<String, Double> sortedLodScores = h.sortHashMapByValuesD(lodScores);

		System.out.println("Hits, sorted low to high:");
		for (String sequenceFeature : sortedLodScores.keySet())
		{
			System.out.println(sequenceFeature + "\t" + sortedLodScores.get(sequenceFeature));
		}

		System.out.println("Writing output dataframe to " + outputFile.getAbsolutePath());
		PrintWriter pw = new PrintWriter(outputFile, "UTF-8");
		pw.println("SNP\tCHR\tBP\tP");
		for (String sequenceFeature : pvals.keySet())
		{
			String loc = h.changeChromLetterToNumber((sequenceFeatureLocations.get(sequenceFeature)));
			pw.println(sequenceFeature + "\t" + loc + "\t" + pvals.get(sequenceFeature));
		}
		pw.flush();
		pw.close();

		// SNP CHR BP P
		// rs1 1 100000 0.9148
		// rs2 1 200000 0.9371
		// rs3 1 400000 0.2861
		// rs4 1 300000 0.8304
		// rs5 1 500000 0.6417
		// rs6 1 600000 0.5191

		String scriptLocation = outputFile.getAbsolutePath() + ".R";
		String plotLocation = outputFile.getAbsolutePath() + ".png";
		System.out.println("Creating plot at " + plotLocation + " using " + scriptLocation);
		h.Rscript(scriptLocation, plotLocation, outputFile.getCanonicalPath());

	}
}
