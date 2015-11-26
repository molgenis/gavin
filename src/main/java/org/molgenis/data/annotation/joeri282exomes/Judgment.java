package org.molgenis.data.annotation.joeri282exomes;

public class Judgment
{
	public enum Classification{
		Benign, Pathogn
		//VOUS,
	}
	
	public enum Confidence{
		High, Medium, Low, FP_FN_25perc
	}
	
	String reason;
	Classification classification;
	Confidence confidence;

	public Judgment(Classification classification, Confidence confidence, String reason)
	{
		super();
		this.reason = reason;
		this.classification = classification;
		this.confidence = confidence;
	}

	public String getReason()
	{
		return reason;
	}
	
	public Classification getClassification()
	{
		return classification;
	}
	
	public Confidence getConfidence()
	{
		return confidence;
	}

	@Override
	public String toString()
	{
		return "Judgment [reason=" + reason + ", classification=" + classification + "]";
	}
	
	
}
