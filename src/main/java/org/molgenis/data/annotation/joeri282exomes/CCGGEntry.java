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
	double PathoMAFThreshold;
	double PopImpactHighPerc;
	double PopImpactModeratePerc;
	double PopImpactLowPerc;
	double PopImpactModifierPerc;
	double PathoImpactHighPerc;
	double PathoImpactModeratePerc;
	double PathoImpactLowPerc;
	double PathoImpactModifierPerc;
	double PopImpactHighEq;
	double PopImpactModerateEq;
	double PopImpactLowEq;
	double PopImpactModifierEq;
	int NrOfCADDScoredPopulationVars;
	int NrOfCADDScoredPathogenicVars;
	double MeanPopulationCADDScore;
	double MeanPathogenicCADDScore;
	double MeanDifference;
	double UTestPvalue;
	double Sens95thPerCADDThreshold;
	double Spec95thPerCADDThreshold;
	String Recommendation;
	
	private enum Category{
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
	}
		
}
