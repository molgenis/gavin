package org.molgenis.data.annotation.joeri282exomes;

import java.util.HashMap;
import java.util.List;

import org.molgenis.data.Entity;

public class CandidateVariant
{
	//record is the base information
	Entity vcfRecord;
	
	//allele and gene further specify the context in which this variant was judged
	String allele;
	int altIndex;
	String gene;
	Double cadd;
	
	//the final classification and reason
	Judgment judgment;
	
	//the samples with interesting genotypes for this variant (+allele/gene)
	HashMap<String, Entity> sampleIds;

	public CandidateVariant(Entity vcfRecord, String allele, int altIndex, String gene, Double cadd, Judgment judgment, HashMap<String, Entity> sampleIds)
	{
		super();
		this.vcfRecord = vcfRecord;
		this.allele = allele;
		this.altIndex = altIndex;
		this.gene = gene;
		this.cadd = cadd;
		this.judgment = judgment;
		this.sampleIds = sampleIds;
	}
	
	public Double getCadd()
	{
		return cadd;
	}

	public int getAltIndex()
	{
		return altIndex;
	}

	public Entity getVcfRecord()
	{
		return vcfRecord;
	}

	public String getAllele()
	{
		return allele;
	}

	public String getGene()
	{
		return gene;
	}

	public Judgment getJudgment()
	{
		return judgment;
	}

	public HashMap<String, Entity> getSampleIds()
	{
		return sampleIds;
	}

	@Override
	public String toString()
	{
		return "CandidateVariant [allele=" + allele + ", gene=" + gene + ", judgment=" + judgment + ", sampleIds="
				+ sampleIds + "]";
	}
	
	
}
