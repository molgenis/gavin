package org.molgenis.caddtlmapping.binom;

import java.io.File;
import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Date;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;

import org.apache.commons.math3.stat.inference.AlternativeHypothesis;
import org.apache.commons.math3.stat.inference.BinomialTest;
import org.molgenis.caddtlmapping.binom.structs.Bin;
import org.molgenis.caddtlmapping.binom.structs.BinnedGenotypeCount;

/**
 * Fresh implementation based on binomial test of over represented genotypes for variants within a certain range of CADD
 * scores. Does not need MAF cutoff to find candidates because we want to use the underlying genotype data to do the
 * statistical test. Takes a folder of CADD files, so we can download multiple and check them all because we want to
 * assess as many variants as possible.
 * 
 * TODO: find out clever way to also use 1000G, GoNL as reference source to enable WGS and not just WES TODO: hemizygous
 * & X/Y PAR ? TODO: multi allelic ok?
 * 
 * how to?
 * 15	66641732	G	C,A,T	AC_Het=11380,1,2,2,0,0;AC_Hom=1094,0,0;
 * 
 * @author jvelde
 *
 */
public class CADDTLMapping
{

	private File vcfFile;
	private File exacFile;
	private File caddFolder;
	private File outputFile;
	private String inheritance;
	private File patientSampleIdsFile;

	public CADDTLMapping(File vcfFile, File exacFile, File caddFolder, String inheritance, File patientSampleIdsFile)
	{
		super();
		this.vcfFile = vcfFile;
		this.exacFile = exacFile;
		this.caddFolder = caddFolder;
		DateFormat dateFormat = new SimpleDateFormat("yyyy-MM-dd-HH-mm-ss");
		Date date = new Date();
		this.outputFile = new File("caddtlmappingresults_" + dateFormat.format(date) + "_" + inheritance + "_" + ".tsv");
		this.inheritance = inheritance;
		this.patientSampleIdsFile = patientSampleIdsFile;
	}

	public void start() throws Exception
	{
		Helper h = new Helper(vcfFile, exacFile, caddFolder);

		// list of sample identifers within the VCF file that represent the patients
		// TODO: immediatly check patient VCF file and give warnings when there are identifiers not mapping etc?
		ArrayList<String> patientSampleIdList = null;
		if (patientSampleIdsFile != null)
		{
			patientSampleIdList = h.loadPatientSampleIdList(patientSampleIdsFile);
			System.out.println("Loaded patient sample identifier list from " + patientSampleIdsFile.getName()
					+ ", found " + patientSampleIdList.size() + " identifiers");
		}

		// scan through ExaC and count eligible genotypes for every sequence feature annoted in "CSQ=", multiple values
		// separated by ",", column 15,
		Bin[] bins = new Bin[]{new Bin(0, 10), new Bin(10, 20), new Bin(20, 30), new Bin(30, 100)};

		HashMap<String, List<BinnedGenotypeCount>> binnedExACGTC = h.getBinnedGenotypeCountsPerSequenceFeature(bins,
				inheritance);

		// scan through patient and count eligible genotypes for every sequence feature annoted in "ANN="
		// TODO: ExAC uses VEP and patient VCF uses SnpEff to get gene names, this is not optimal, might be
		// differences.. could use SnpEff on ExAC perhaps? or VEP on the patients? or just report the mismatches?
		HashMap<String, List<BinnedGenotypeCount>> binnedPatientGTC = h.getBinnedPatientGenotypeCountsPerSequenceFeature(bins, inheritance, patientSampleIdList);

		
		HashMap<String, Double> lodScores = new HashMap<String, Double>();
		
		// now perform binomial test per sequence feature, per bin (so essentially multiple tests per gene, but not
		// many)
		for (String sequenceFeature : binnedPatientGTC.keySet())
		{
			if (!binnedExACGTC.keySet().contains(sequenceFeature))
			{
				System.out.println("WARNING: ExAC does not have a gene named '" + sequenceFeature + "' ! skipping..");
				continue;
			}

			for (BinnedGenotypeCount patientGTC : binnedPatientGTC.get(sequenceFeature))
			{
				System.out.println(sequenceFeature + ", cadd " + patientGTC.bin.lower + "-" + patientGTC.bin.upper);
				int successPat = patientGTC.actingGenotypes;
				int drawsPat = patientGTC.totalGenotypes;

				double prob = -1;
				for (BinnedGenotypeCount exacGTC : binnedPatientGTC.get(sequenceFeature))
				{
					if (exacGTC.bin.lower == patientGTC.bin.lower && exacGTC.bin.upper == patientGTC.bin.upper)
					{
						double successExAC = patientGTC.actingGenotypes;
						double drawsExAC = patientGTC.totalGenotypes;
						prob = successExAC / drawsExAC;
						
					}
				}

				BinomialTest binom = new BinomialTest();
				double pval = binom.binomialTest(drawsPat, successPat, prob, AlternativeHypothesis.GREATER_THAN);
				double lod = -Math.log10(pval);
		//		System.out.println(sequenceFeature + " has LOD " + lod);
				lodScores.put(sequenceFeature, lod);
			}

		}
		
		
		
		System.out.println("Processing the results..");
		LinkedHashMap<String, Double> sortedLodScores = h.sortHashMapByValuesD(lodScores);

		System.out.println("Hits, sorted low to high:");
		for (String sequenceFeature : sortedLodScores.keySet())
		{
			System.out.println(sequenceFeature + "\t" + sortedLodScores.get(sequenceFeature));
		}
		
		

	}
}
