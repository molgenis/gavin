package org.molgenis.data.annotation.joeri282exomes.struct;

/**
 * 
 * We store counts for pathogenic and VOUS variants in 3 values each : total allele count, recessive genotype count, dominant genotype count.
 * 
 * toString examples:
 * 
 * 2,0,2|23,5,13
 * 5,1,3|42,10,22
 * 
 * 
 * @author jvelde
 *
 */
public class AlleleCounts
{

	public static String printNull = "0,0,0|0,0,0";
	
	public int pathoTotalAlleles;
	public int pathoRecessiveGenotypes;
	public int pathoDominantGenotypes;
	public int vousTotalAlleles;
	public int vousRecessiveGenotypes;
	public int vousDominantGenotypes;
	
	public AlleleCounts()
	{
		super();
		this.pathoTotalAlleles = 0;
		this.pathoRecessiveGenotypes = 0;
		this.pathoDominantGenotypes = 0;
		this.vousTotalAlleles = 0;
		this.vousRecessiveGenotypes = 0;
		this.vousDominantGenotypes = 0;
		
	}

	@Override
	public String toString()
	{
		return pathoTotalAlleles+","+pathoRecessiveGenotypes+","+pathoDominantGenotypes+"|"+vousTotalAlleles+","+vousRecessiveGenotypes+","+vousDominantGenotypes;
	}
	
	

}
