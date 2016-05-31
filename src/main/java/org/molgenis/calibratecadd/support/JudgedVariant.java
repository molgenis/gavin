package org.molgenis.calibratecadd.support;

import org.molgenis.data.Entity;
import org.molgenis.data.annotation.entity.impl.gavin.Judgment;


public class JudgedVariant
{
	Judgment judgment;
	Entity e;
	
	public enum ExpertClassification{
		B, LB, V, LP, P
	}
	
	ExpertClassification expertClassification;
	
	public Judgment getJudgment()
	{
		return judgment;
	}
	public Entity getE()
	{
		return e;
	}
	public ExpertClassification getExpertClassification()
	{
		return expertClassification;
	}
	public JudgedVariant(Judgment judgment, Entity variant, ExpertClassification expertClassification)
	{
		super();
		this.judgment = judgment;
		this.e = variant;
		this.expertClassification = expertClassification;
	}
	
	public String printVariant()
	{
		return e.getString("ID") + " (" + e.getString("#CHROM") + ":" + e.getString("POS") + " " + e.getString("REF") + " to " + e.getString("ALT") + "), reason: " + judgment.getReason();
	}
	@Override
	public String toString()
	{
		return "JudgedVariant [judgment=" + judgment + ", e=" + e + ", expertClassification=" + expertClassification
				+ ", getJudgment()=" + getJudgment() + ", getE()=" + getE() + ", getExpertClassification()="
				+ getExpertClassification() + ", printVariant()=" + printVariant() + ", getClass()=" + getClass()
				+ ", hashCode()=" + hashCode() + ", toString()=" + super.toString() + "]";
	}
	
	
	
}
