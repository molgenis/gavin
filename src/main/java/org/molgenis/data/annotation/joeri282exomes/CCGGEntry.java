package org.molgenis.data.annotation.joeri282exomes;

public class CCGGEntry
{
	String gene;
	Category category;
	String chromosome;
	int start;
	int end;
	int NrOfPopulationVariants;
	int NrOfPathogenicVariants;
	int NrOfOverlappingVariants;
	int NrOfFilteredPopVariants;
	Double PathoMAFThreshold;
	Double PopImpactHighPerc;
	Double PopImpactModeratePerc;
	Double PopImpactLowPerc;
	Double PopImpactModifierPerc;
	Double PathoImpactHighPerc;
	Double PathoImpactModeratePerc;
	Double PathoImpactLowPerc;
	Double PathoImpactModifierPerc;
	Double PopImpactHighEq;
	Double PopImpactModerateEq;
	Double PopImpactLowEq;
	Integer PopImpactModifierEq;
	Integer NrOfCADDScoredPopulationVars;
	Integer NrOfCADDScoredPathogenicVars;
	Double MeanPopulationCADDScore;
	Double MeanPathogenicCADDScore;
	Double MeanDifference;
	Double UTestPvalue;
	Double Sens95thPerCADDThreshold;
	Double Spec95thPerCADDThreshold;
	String Recommendation;
	
	enum Category{
		N1, N2, T1, T2, I1, I2, I3, C1, C2, C3, C4, C5
	}
	
	public CCGGEntry(String lineFromFile) throws Exception
	{
		String[] split = lineFromFile.split("\t", -1);
		if(split.length != 31)
		{
			throw new Exception("not 31 elements");
		}
		
		this.gene = split[0];
		this.category = Category.valueOf(split[1]);
		this.chromosome = split[2];
		this.start = Integer.valueOf(split[3]);
		this.end = Integer.valueOf(split[4]);
		this.NrOfPopulationVariants = Integer.valueOf(split[5]);
		this.NrOfPathogenicVariants = Integer.valueOf(split[6]);
		this.NrOfOverlappingVariants = Integer.valueOf(split[7]);
		this.NrOfFilteredPopVariants = Integer.valueOf(split[8]);
		this.PathoMAFThreshold = split[9].isEmpty() ? null : Double.parseDouble(split[9]);
		this.MeanPopulationCADDScore = split[24].isEmpty() ? null : Double.parseDouble(split[24]);
		this.MeanPathogenicCADDScore = split[25].isEmpty() ? null : Double.parseDouble(split[25]);
		this.MeanDifference = split[26].isEmpty() ? null : Double.parseDouble(split[26]);
		this.UTestPvalue = split[27].isEmpty() ? null : Double.parseDouble(split[27]);
		this.Sens95thPerCADDThreshold = split[28].isEmpty() ? null : Double.parseDouble(split[28]);
		this.Spec95thPerCADDThreshold = split[29].isEmpty() ? null : Double.parseDouble(split[29]);
	}
		
}
