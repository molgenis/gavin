package org.molgenis.data.annotation.joeri282exomes;

public class Judgment
{
	public enum Classification{
		Benign, Likely_Benign, Likely_Pathogenic, Pathogenic
		//VOUS,
	}
	
	String reason;
	Classification classification;

	public Judgment(Classification classification, String reason)
	{
		super();
		this.reason = reason;
		this.classification = classification;
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
