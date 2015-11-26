package org.molgenis.calibratecadd.support;

import org.molgenis.data.Entity;
import org.molgenis.data.annotation.joeri282exomes.Judgment;

public class JudgedVariant
{
	Judgment j;
	Entity e;
	
	public enum ExpertClassification{
		B, LB, V, LP, P
	}
	
	ExpertClassification expertClassification;
	
	public Judgment getJudgment()
	{
		return j;
	}
	public Entity getE()
	{
		return e;
	}
	public ExpertClassification getExpertClassification()
	{
		return expertClassification;
	}
	public JudgedVariant(Judgment j, Entity e, ExpertClassification expertClassification)
	{
		super();
		this.j = j;
		this.e = e;
		this.expertClassification = expertClassification;
	}
	
	
}
