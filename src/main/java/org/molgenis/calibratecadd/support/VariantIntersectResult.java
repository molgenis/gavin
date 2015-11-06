package org.molgenis.calibratecadd.support;

import java.util.List;

public class VariantIntersectResult
{
	public List<EntityPlus> inBoth_exac;
	public List<EntityPlus> inBoth_clinvar;
	public List<EntityPlus> inClinVarOnly;
	public List<EntityPlus> inExACOnly;
	
	
	public VariantIntersectResult(List<EntityPlus> inBoth_exac, List<EntityPlus> inBoth_clinvar,
			List<EntityPlus> inClinVarOnly, List<EntityPlus> inExACOnly)
	{
		super();
		this.inBoth_exac = inBoth_exac;
		this.inBoth_clinvar = inBoth_clinvar;
		this.inClinVarOnly = inClinVarOnly;
		this.inExACOnly = inExACOnly;
	}

	
	
	
	
}
