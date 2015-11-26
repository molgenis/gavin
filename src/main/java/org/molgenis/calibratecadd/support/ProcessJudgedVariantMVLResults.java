package org.molgenis.calibratecadd.support;

import java.util.HashMap;
import java.util.List;

import org.molgenis.calibratecadd.support.JudgedVariant.ExpertClassification;
import org.molgenis.data.annotation.joeri282exomes.Judgment.Classification;
import org.molgenis.data.annotation.joeri282exomes.Judgment.Confidence;

public class ProcessJudgedVariantMVLResults
{
	public static void printResults(HashMap<String, List<JudgedVariant>> judgedMVLVariants)
	{
//		printCountsOfExpertMVLClassifications(judgedMVLVariants);
//		printCountsOfCCGGMVLClassifications(judgedMVLVariants, Confidence.FP_FN_1perc);
//		printCountsOfCCGGMVLClassifications(judgedMVLVariants, Confidence.FP_FN_5perc);
		calculateAndPrint_FP_FN_stats(judgedMVLVariants, Confidence.Medium);
	}
	
	public static void calculateAndPrint_FP_FN_stats(HashMap<String, List<JudgedVariant>> judgedMVLVariants, Confidence confidenceTranche)
	{
		System.out.println("\nFalse posities & false negatives in confidence tranche: " + confidenceTranche);
		int grandTotal = 0;
		
		System.out.println("\t" + "TN" + "\t" + "TP" + "\t" + "FP" + "\t" + "FN");
		
		for(String mvl : judgedMVLVariants.keySet())
		{
			int TN = 0;
			int TP = 0;
			int FP = 0;
			int FN = 0;
			for(JudgedVariant jv : judgedMVLVariants.get(mvl))
			{
				if(jv.getJudgment() != null && jv.getJudgment().getConfidence().equals(confidenceTranche))
				{
					if((jv.getExpertClassification().equals(ExpertClassification.B) || jv.getExpertClassification().equals(ExpertClassification.LB)) && jv.getJudgment().getClassification().equals(Classification.Benign))
					{
						TN ++;
					}
					else if((jv.getExpertClassification().equals(ExpertClassification.P) || jv.getExpertClassification().equals(ExpertClassification.LP)) && jv.getJudgment().getClassification().equals(Classification.Pathogn))
					{
						TP ++;
					}
					else if((jv.getExpertClassification().equals(ExpertClassification.B) || jv.getExpertClassification().equals(ExpertClassification.LB)) && jv.getJudgment().getClassification().equals(Classification.Pathogn))
					{
						FP ++;
					}
					else if((jv.getExpertClassification().equals(ExpertClassification.P) || jv.getExpertClassification().equals(ExpertClassification.LP)) && jv.getJudgment().getClassification().equals(Classification.Benign))
					{
						FN ++;
					}
				}
			}
			
			System.out.println(mvl + "\t" + TN + "\t" + TP + "\t" + FP + "\t" + FN);
			
			
		}
		

	}

	public static void printCountsOfCCGGMVLClassifications(HashMap<String, List<JudgedVariant>> judgedMVLVariants, Confidence confidenceTranche)
	{
		System.out.println("\nCounts in confidence tranche: " + confidenceTranche);
		int grandTotal = 0;
		for(String mvl : judgedMVLVariants.keySet())
		{
			System.out.print("\t" + mvl);
		}
		System.out.println("\t" + "TOTAL");
		int totalPerCL = 0;
		for(Classification cl : Classification.values())
		{
			System.out.print(cl + "\t");
			int count = 0;
			for(String mvl : judgedMVLVariants.keySet())
			{	
				for(JudgedVariant jv : judgedMVLVariants.get(mvl))
				{
					if(jv.getJudgment() != null && jv.getJudgment().getClassification().equals(cl) && jv.getJudgment().getConfidence().equals(confidenceTranche))
					{
						count++;
					}
				}
				System.out.print(count + "\t");
				totalPerCL += count;
				count = 0;
			}
			grandTotal += totalPerCL;
			System.out.println(totalPerCL);
			totalPerCL = 0;
		}
		System.out.print("TOTAL");
		int count = 0;
		for(String mvl : judgedMVLVariants.keySet())
		{	
			for(JudgedVariant jv : judgedMVLVariants.get(mvl))
			{
				if(jv.getJudgment() != null && jv.getJudgment().getConfidence().equals(confidenceTranche))
				{
					count++;
				}
			}
			System.out.print("\t" + count);
			count = 0;
		}
		System.out.println("\t" + grandTotal);
	}
	
	

	public static void printCountsOfExpertMVLClassifications(HashMap<String, List<JudgedVariant>> judgedMVLVariants)
	{
		System.out.println("\nExpert MVL classifications");
		int grandTotal = 0;
		for(String mvl : judgedMVLVariants.keySet())
		{
			System.out.print("\t" + mvl);
		}
		System.out.println("\t" + "TOTAL");
		int totalPerEC = 0;
		for(ExpertClassification ec : ExpertClassification.values())
		{
			System.out.print(ec + "\t");
			int count = 0;
			for(String mvl : judgedMVLVariants.keySet())
			{	
				for(JudgedVariant jv : judgedMVLVariants.get(mvl))
				{
					if(jv.getExpertClassification().equals(ec))
					{
						count++;
					}
				}
				System.out.print(count + "\t");
				totalPerEC += count;
				count = 0;
			}
			grandTotal += totalPerEC;
			System.out.println(totalPerEC);
			totalPerEC = 0;
		}
		System.out.print("TOTAL");
		int count = 0;
		for(String mvl : judgedMVLVariants.keySet())
		{	
			for(JudgedVariant jv : judgedMVLVariants.get(mvl))
			{
				count++;
			}
			System.out.print("\t" + count);
			count = 0;
		}
		System.out.println("\t" + grandTotal);
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
