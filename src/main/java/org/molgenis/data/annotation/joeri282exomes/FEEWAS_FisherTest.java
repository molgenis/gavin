package org.molgenis.data.annotation.joeri282exomes;

import java.io.File;
import java.util.HashMap;

import org.molgenis.data.annotation.joeri282exomes.legacy.DifferentialExomes;
import org.molgenis.data.annotation.joeri282exomes.legacy.FishersExactTest;
import org.molgenis.data.annotation.joeri282exomes.struct.AlleleCounts;
import org.molgenis.data.annotation.joeri282exomes.struct.GeneGroupsAlleleCountUtils;

public class FEEWAS_FisherTest
{
	
	public FEEWAS_FisherTest(File feewasFile) throws Exception
	{
		GeneGroupsAlleleCountUtils data = new GeneGroupsAlleleCountUtils();
		data.readFromFile(feewasFile);
		
//		System.out.println("data read");
		
		//gene to group, group to pval
		HashMap<String, HashMap<String, Double>> pvalsDom = new HashMap<String, HashMap<String, Double>>();
		HashMap<String, HashMap<String, Double>> pvalsAdd = new HashMap<String, HashMap<String, Double>>();
		HashMap<String, HashMap<String, Double>> pvalsRec = new HashMap<String, HashMap<String, Double>>();
		
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
				
				HashMap<String, Double> genePvalsDom = pvalsDom.get(gene);
				if(genePvalsDom == null) {
					genePvalsDom = new HashMap<String, Double>();
					pvalsDom.put(gene, genePvalsDom);
				}
				genePvalsDom.put(groupToBeTested, pvalDominant);
				
				HashMap<String, Double> genePvalsAdd = pvalsAdd.get(gene);
				if(genePvalsAdd == null) {
					genePvalsAdd = new HashMap<String, Double>();
					pvalsAdd.put(gene, genePvalsAdd);
				}
				genePvalsAdd.put(groupToBeTested, pvalAdditive);
				
				HashMap<String, Double> genePvalsRec = pvalsRec.get(gene);
				if(genePvalsRec == null) {
					genePvalsRec = new HashMap<String, Double>();
					pvalsRec.put(gene, genePvalsRec);
				}
				genePvalsRec.put(groupToBeTested, pvalRecessive);

				
				double lodDominant = -Math.log10(pvalDominant);
				double lodAdditive = -Math.log10(pvalAdditive);
				double lodRecessive = -Math.log10(pvalDominant);
				
				if(lodDominant > 5)
				{
					System.out.println("LOD score for dominant model " + lodDominant + " for group " + groupToBeTested + " in gene " + gene);
				//	System.out.println("inputs for fisher test: " + ac.getActingDominant() + ", " + ac.getNonActingDominant() + ", " + allOtherGroups.getActingDominant() + ", " + allOtherGroups.getNonActingDominant());
				}
				
				if(lodAdditive > 5)
				{
					System.out.println("LOD score for additive model " + lodAdditive + " for group " + groupToBeTested + " in gene " + gene);
				}
				
				if(lodRecessive > 5)
				{
					System.out.println("LOD score for recessive model " + lodRecessive + " for group " + groupToBeTested + " in gene " + gene);
				}
				
			}
		}
		
		File outputDataFrameFileDom = new File(feewasFile.getAbsolutePath() + "_dataframe_dom");
		File outputDataFrameFileAdd = new File(feewasFile.getAbsolutePath() + "_dataframe_add");
		File outputDataFrameFileRec = new File(feewasFile.getAbsolutePath() + "_dataframe_rec");
		DifferentialExomes.writeDataSetForPlotScript(outputDataFrameFileDom, data.getGeneLocs(), data.getGroups(), pvalsDom);
		DifferentialExomes.writeDataSetForPlotScript(outputDataFrameFileAdd, data.getGeneLocs(), data.getGroups(), pvalsAdd);
		DifferentialExomes.writeDataSetForPlotScript(outputDataFrameFileRec, data.getGeneLocs(), data.getGroups(), pvalsRec);
		DifferentialExomes.plot(outputDataFrameFileDom, data.getGroups());
		DifferentialExomes.plot(outputDataFrameFileAdd, data.getGroups());
		DifferentialExomes.plot(outputDataFrameFileRec, data.getGroups());
		
	}

	public static void main(String[] args) throws Exception
	{
		if(args.length != 1)
		{
			throw new Exception("need to supply 1 arg: FEEWAS file");
		}
		File feewasFile = new File(args[0]);
		new FEEWAS_FisherTest(feewasFile);

	}

}
