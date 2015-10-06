package org.molgenis.calibratecadd.structs;

import java.util.List;

import org.molgenis.data.Entity;

public class VariantIntersectResult
{
	public List<Entity> inBoth_exac;
	public List<Entity> inBoth_clinvar;
	public List<Entity> inClinVarOnly;
	public List<Entity> inExACOnly;
	
	public VariantIntersectResult(List<Entity> inBoth_exac, List<Entity> inBoth_clinvar, List<Entity> inClinVarOnly,
			List<Entity> inExACOnly)
	{
		super();
		this.inBoth_exac = inBoth_exac;
		this.inBoth_clinvar = inBoth_clinvar;
		this.inClinVarOnly = inClinVarOnly;
		this.inExACOnly = inExACOnly;
	}
	
	
	
	
	
}
