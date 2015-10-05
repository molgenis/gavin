package org.molgenis.calibratecadd.structs;

import java.util.List;

import org.molgenis.data.Entity;

public class VariantIntersectResult
{
	public List<Entity> inBoth;
	public List<Entity> inClinVarOnly;
	public List<Entity> inExACOnly;
	
	
	public VariantIntersectResult(List<Entity> inBoth, List<Entity> inClinVarOnly, List<Entity> inExACOnly)
	{
		super();
		this.inBoth = inBoth;
		this.inClinVarOnly = inClinVarOnly;
		this.inExACOnly = inExACOnly;
	}
	
	
}
