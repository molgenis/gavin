package org.molgenis.caddtlmapping;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;

import org.molgenis.data.Entity;
import org.molgenis.data.vcf.VcfRepository;

// "CADD-TL mapping"
// test 3 inheritance models: recessive, additive, dominant
// test 3 frequency models: singleton (0%), rare (<1%), common (<5%)
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
	private File caddFile;
	File outputFile;
	private double mafThreshold;
	private String inheritance;
	private File patientSampleIdsFile;

	public CADDTLMapping(File vcfFile, File exacFile, File caddFile, File outputFile, double mafThreshold,
			String inheritance, File patientSampleIdsFile)
	{
		super();
		this.vcfFile = vcfFile;
		this.exacFile = exacFile;
		this.caddFile = caddFile;
		this.outputFile = outputFile;
		this.mafThreshold = mafThreshold;
		this.inheritance = inheritance;
		this.patientSampleIdsFile = patientSampleIdsFile;
	}

	/**
	 * Central function of the tool
	 * 
	 * @throws IOException
	 */
	public void start() throws IOException
	{

		// list of sample identifers within the VCF file that represent the patients
		// TODO: check VCF file and give warnings when there are identifiers not mapping etc
		System.out.println("Loading patient sampleID list..");
		ArrayList patientSampleIdList = Helper.loadPatientSampleIdList(patientSampleIdsFile);

		System.out.println("First pass to identify candidate variants per sequence feature..");
		HashMap<String, List<Integer>> candidates = Helper.getCandidateVariantsPerGene(vcfFile, exacFile, mafThreshold);

	}
}
