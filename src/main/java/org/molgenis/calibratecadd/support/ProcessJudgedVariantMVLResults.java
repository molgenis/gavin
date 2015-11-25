package org.molgenis.calibratecadd.support;

import java.util.HashMap;
import java.util.List;

public class ProcessJudgedVariantMVLResults
{
	public static void printResults(HashMap<String, List<JudgedVariant>> judgedMVLVariants)
	{
		System.out.println("MVL:");
		for(String mvl : judgedMVLVariants.keySet())
		{
			System.out.println(mvl);
		}
	}
}

//System.out.println("\n\nTotal amount of B/LB/LP/P variants in original MVLs = " + nrOfMVL_B_LB_P_LP_variants + " ("+nrOfMVL_B_LB_variants+" benign, "+nrOfMVL_P_LP_variants+" pathogenic)");
//System.out.println("Total amount of B/LB/LP/P variants that CCGG-classified as B/LB/LP/P = " + nrOf_B_LB_P_LP_CCGGclassifiedVariants_acrossAllMVLs);
//
//System.out.println("\nCorrectly classified benign = " + correctlyClassifiedBenign);
//System.out.println("Wrongly classified benign (LB/B judged as LP/P, false positives) = " + wronglyClassifiedBenign);
//System.out.println("Correctly classified pathogenic = " + correctlyClassifiedPathogenic);
//System.out.println("Wrongly classified pathogenic (LP/P judged as LB/B, false negatives) = " + wronglyClassifiedPathogenic);
//
//System.out.println("\nAverage false positive rate in CCGG-classifications across all MVLs = " + Math.round(((double)wronglyClassifiedBenign/(double)(wronglyClassifiedBenign+correctlyClassifiedBenign))*100) + "%");
//System.out.println("Average false negative rate in CCGG-classifications across all MVLs = " + Math.round(((double)wronglyClassifiedPathogenic/(double)(wronglyClassifiedPathogenic+correctlyClassifiedPathogenic))*100) + "%");
//
//System.out.println("\nClassification of VOUS variants:");
//System.out.println("Total amount of VOUS variants in original MVLs = " + nrOfMVL_VOUS_variants);
//System.out.println("Total amount of VOUS variants that CCGG-classified as B/LB/LP/P = " + nrOf_VOUS_CCGGclassifiedVariants_acrossAllMVLs + " (" + nrOf_VOUS_B_LB_CCGGclassifiedVariants_acrossAllMVLs + " benign, " + nrOf_VOUS_P_LP_CCGGclassifiedVariants_acrossAllMVLs + " pathogenic)");
