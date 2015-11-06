package org.molgenis.calibratecadd.support;

public class ClinVarVariant
{

	public String chrom;
	public long pos;
	public String id;
	public String ref;
	public String alt;
	
	public String name;
	public String gene;
	public String clinsig;
	
	public String clinvarInfoToString()
	{
		return name + "|" + gene + "|" + clinsig.replace(";", "/");
	}
	
	
	public ClinVarVariant(String chrom, String pos, String id, String ref, String alt, String name, String gene,
			String clinsig)
	{
		super();
		this.chrom = chrom;
		this.pos = Long.parseLong(pos);
		this.id = id;
		this.ref = ref;
		this.alt = alt;
		this.name = name;
		this.gene = gene;
		this.clinsig = clinsig;
	}
	
	
	
}
