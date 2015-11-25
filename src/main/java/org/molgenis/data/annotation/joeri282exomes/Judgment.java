package org.molgenis.data.annotation.joeri282exomes;

public class Judgment
{
	public enum Classification{
		Benign, Pathogenic
		//VOUS,
	}
	
	public enum Confidence{
		FP_FN_1perc, FP_FN_5perc, FP_FN_10perc, FP_FN_25perc
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

	@Override
	public String toString()
	{
		return "Judgment [reason=" + reason + ", classification=" + classification + "]";
	}
	
	
}
