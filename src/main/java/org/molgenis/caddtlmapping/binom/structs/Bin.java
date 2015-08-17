package org.molgenis.caddtlmapping.binom.structs;

public class Bin
{
	public Bin(double lower, double upper)
	{
		this.lower = lower;
		this.upper = upper;
	}
	public double lower;
	public double upper;
	
	@Override
	public String toString()
	{
		return "Bin [lower=" + lower + ", upper=" + upper + "]";
	}
	
	
}
