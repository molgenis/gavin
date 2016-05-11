package org.molgenis.calibratecadd.support;

public class BootStrappingVariant
{
	public enum OutCome{
		TP, TN, FP, FN, VOUS
	}

	public enum ExpClsf{
		B, P
	}

	public enum Label{
		C1_C2, C4, C3
	}
	
	private Label label;
	private OutCome outcome;
	private ExpClsf expClsf;

	public BootStrappingVariant(OutCome outcome, Label label, ExpClsf expClsf)
	{
		super();
		this.outcome = outcome;
		this.label = label;
		this.expClsf = expClsf;
	}

	public Label getLabel() {
		return label;
	}

	public OutCome getOutcome()
	{
		return outcome;
	}

	public ExpClsf getExpClsf() { return expClsf; }

	@Override
	public String toString() {
		return "BootStrappingVariant{" +
				"label=" + label +
				", outcome=" + outcome +
				", expClsf=" + expClsf +
				'}';
	}
}
