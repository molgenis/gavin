package org.molgenis.data.annotation.joeri282exomes;

public class Judgment
{
	public enum Classification{
		Benign, Pathogn, VOUS
	}
	
	public enum Method{
		calibrated, genomewide
	}
	
	String reason;
	Classification classification;
	Method method;
	String gene;

	public Judgment(Classification classification, Method method, String gene, String reason)
	{
		super();
		this.reason = reason;
		this.gene = gene;
		this.classification = classification;
		this.method = method;
	}

	public String getReason()
	{
		return reason;
	}
	
	public Classification getClassification()
	{
		return classification;
	}
	
	public Method getConfidence()
	{
		return method;
	}

	@Override
	public String toString()
	{
		return "Judgment [reason=" + reason + ", classification=" + classification + "]";
	}

	public String getGene() {
		return gene;
	}
}
