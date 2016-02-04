package org.molgenis.data.annotation.joeri282exomes;

import java.io.FileNotFoundException;

import org.molgenis.data.annotation.joeri282exomes.legacy.FishersExactTest;
import org.molgenis.data.annotation.joeri282exomes.struct.AlleleCounts;
import org.molgenis.data.annotation.joeri282exomes.struct.GeneGroupsAlleleCountUtils;

public class FEEWAS_FisherTest
{
	
	public FEEWAS_FisherTest() throws Exception
	{
		GeneGroupsAlleleCountUtils data = new GeneGroupsAlleleCountUtils();
		data.readFromFile("/Users/jvelde/FEEWAS_Patho.tsv");
		
//		System.out.println("data read");
		
		for(String gene : data.getGenesToGroupsToAlleleCounts().keySet())
		{
//			System.out.println("gene " + gene);
			for(String groupToBeTested : data.getGroups())
			{
				AlleleCounts ac = data.getGenesToGroupsToAlleleCounts().get(gene).get(groupToBeTested);
				
				AlleleCounts allOtherGroups = new AlleleCounts();
				for(String otherGroup : data.getGroups())
				{
					if(!otherGroup.equals(groupToBeTested))
					{
						//add up data from all other groups
						AlleleCounts otherGroupAC = data.getGenesToGroupsToAlleleCounts().get(gene).get(otherGroup);
						allOtherGroups.add(otherGroupAC);
					}
				}
				
//				System.out.println("testing " + gene + " in group " + groupToBeTested + " vs others");
				
				//test AC in the 3 different models (dominant, additive, recessive)
				//a = acting, group to be tested
				//b = nonacting, group to the tested
				//c = acting, all other groups
				//d = nonacting, all other groups
				double pvalDominant = FishersExactTest.run(ac.getActingDominant(), ac.getNonActingDominant(), allOtherGroups.getActingDominant(), allOtherGroups.getNonActingDominant());
				double pvalAdditive = FishersExactTest.run(ac.getActingAdditive(), ac.getNonActingAdditive(), allOtherGroups.getActingAdditive(), allOtherGroups.getNonActingAdditive());
				double pvalRecessive = FishersExactTest.run(ac.getActingRecessive(), ac.getNonActingRecessive(), allOtherGroups.getActingRecessive(), allOtherGroups.getNonActingRecessive());
					
				double lodDominant = -Math.log10(pvalDominant);
				double lodAdditive = -Math.log10(pvalAdditive);
				double lodRecessive = -Math.log10(pvalDominant);
				
				if(lodDominant > 3)
				{
					System.out.println("LOD score for dominant model " + lodDominant + " for group " + groupToBeTested + " in gene " + gene);
				//	System.out.println("inputs for fisher test: " + ac.getActingDominant() + ", " + ac.getNonActingDominant() + ", " + allOtherGroups.getActingDominant() + ", " + allOtherGroups.getNonActingDominant());
				}
				
				if(lodAdditive > 3)
				{
					System.out.println("LOD score for additive model " + lodAdditive + " for group " + groupToBeTested + " in gene " + gene);
				}
				
				if(lodRecessive > 3)
				{
					System.out.println("LOD score for recessive model " + lodRecessive + " for group " + groupToBeTested + " in gene " + gene);
				}
				
			}
		}
		
		
	}

	public static void main(String[] args) throws Exception
	{
		new FEEWAS_FisherTest();

	}

}
