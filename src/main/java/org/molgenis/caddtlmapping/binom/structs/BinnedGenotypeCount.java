package org.molgenis.caddtlmapping.binom.structs;

public class BinnedGenotypeCount
{
	public BinnedGenotypeCount(Bin bin, int actingGenotypes, int totalGenotypes)
	{
		this.bin = bin;
		this.actingGenotypes = actingGenotypes;
		this.totalGenotypes = totalGenotypes;
	}

	public Bin bin;
	public int actingGenotypes;
	public int totalGenotypes;
	
	@Override
	public String toString()
	{
		return "BinnedGenotypeCount [bin=" + bin + ", actingGenotypes=" + actingGenotypes + ", totalGenotypes="
				+ totalGenotypes + "]";
	}
	
	
}
