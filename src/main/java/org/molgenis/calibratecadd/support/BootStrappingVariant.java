package org.molgenis.calibratecadd.support;

public class BootStrappingVariant
{
	public enum OutCome{
		TP, TN, FP, FN, VOUS
	}
	
	boolean inCalibGene;
	OutCome outcome;

	public BootStrappingVariant(OutCome outcome, boolean inCalibGene)
	{
		super();
		this.outcome = outcome;
		this.inCalibGene = inCalibGene;
	}

	public boolean isInCalibGene()
	{
		return inCalibGene;
	}

	public OutCome getOutcome()
	{
		return outcome;
	}
	
	
}
