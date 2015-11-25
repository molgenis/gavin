package org.molgenis.calibratecadd.support;

import org.molgenis.data.Entity;
import org.molgenis.data.annotation.joeri282exomes.Judgment;

public class JudgedVariant
{
	Judgment j;
	Entity e;
	String expertClassification;
	
	public Judgment getJ()
	{
		return j;
	}
	public Entity getE()
	{
		return e;
	}
	public String getExpertClassification()
	{
		return expertClassification;
	}
	public JudgedVariant(Judgment j, Entity e, String expertClassification)
	{
		super();
		this.j = j;
		this.e = e;
		this.expertClassification = expertClassification;
	}
	
	
}
