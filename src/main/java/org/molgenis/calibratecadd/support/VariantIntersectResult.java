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


	@Override
	public String toString() {

		StringBuffer sb = new StringBuffer();
		sb.append("Variant Intersection Result:\n");
		sb.append("- In both sets, ExAC data:\n");
		for(EntityPlus e : inBoth_exac)
		{
			sb.append("\t"+e.toString() + "\n");
		}

		sb.append("- In both sets, ClinVar data:\n");
		for(EntityPlus e : inBoth_clinvar)
		{
			sb.append("\t"+e.toString() + "\n");
		}

		sb.append("- In ClinVar only:\n");
		for(EntityPlus e : inClinVarOnly)
		{
			sb.append("\t"+e.toString() + "\n");
		}

		sb.append("- In ExAC only:\n");
		for(EntityPlus e : inExACOnly)
		{
			sb.append("\t"+e.toString() + "\n");
		}

		return sb.toString();
	}
}
