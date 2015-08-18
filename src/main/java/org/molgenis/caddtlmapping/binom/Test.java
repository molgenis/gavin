package org.molgenis.caddtlmapping.binom;

import org.apache.commons.math3.stat.inference.AlternativeHypothesis;
import org.apache.commons.math3.stat.inference.BinomialTest;

public class Test
{

	public static void main(String[] args)
	{
		BinomialTest binom = new BinomialTest();
		double pval = binom.binomialTest(201, 5, (19162.0 / 8385655.0), AlternativeHypothesis.GREATER_THAN);
		double lod = -Math.log10(pval);
		System.out.println("pval: " + pval + ", lod: " + lod);
	}
}
