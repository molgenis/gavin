package org.molgenis.data.annotation.joeri282exomes;

import java.util.List;

import org.molgenis.data.Entity;

public class CandidateVariant
{
	//record is the base information
	Entity vcfRecord;
	
	//allele and gene further specify the context in which this variant was judged
	String allele;
	String gene;
	
	//the final classification and reason
	Judgment judgment;
	
	//the samples with interesting genotypes for this variant (+allele/gene)
	List<String> sampleIds;

	public CandidateVariant(Entity vcfRecord, String allele, String gene, Judgment judgment, List<String> sampleIds)
	{
		super();
		this.vcfRecord = vcfRecord;
		this.allele = allele;
		this.gene = gene;
		this.judgment = judgment;
		this.sampleIds = sampleIds;
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

	public List<String> getSampleIds()
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
