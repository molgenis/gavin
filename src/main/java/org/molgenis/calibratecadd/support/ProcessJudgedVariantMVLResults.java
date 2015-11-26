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
		printCountsOfExpertMVLClassifications(judgedMVLVariants);
		printCountsOfCCGGMVLClassifications(judgedMVLVariants, Confidence.High);
		calculateAndPrint_FP_FN_stats(judgedMVLVariants, Confidence.High);
		printCountsOfCCGGMVLClassifications(judgedMVLVariants, Confidence.Medium);
		calculateAndPrint_FP_FN_stats(judgedMVLVariants, Confidence.Medium);
		printCountsOfCCGGMVLClassifications(judgedMVLVariants, Confidence.Low);
		calculateAndPrint_FP_FN_stats(judgedMVLVariants, Confidence.Low);
	}
	
	public static void calculateAndPrint_FP_FN_stats(HashMap<String, List<JudgedVariant>> judgedMVLVariants, Confidence confidenceTranche)
	{
		System.out.println("\nFalse posities & false negatives in confidence tranche: " + confidenceTranche);
		int grandTotal = 0;
		
		System.out.println("\t" + "#TN" + "\t" + "#TP" + "\t" + "#FP" + "\t" + "#FN" + "\t" + "TPR" + "\t" + "TNR" + "\t" + "PPV" + "\t" + "NPV");
		
		int TN_all = 0;
		int TP_all = 0;
		int FP_all = 0;
		int FN_all = 0;
		
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
						TN_all ++;
					}
					else if((jv.getExpertClassification().equals(ExpertClassification.P) || jv.getExpertClassification().equals(ExpertClassification.LP)) && jv.getJudgment().getClassification().equals(Classification.Pathogn))
					{
						TP ++;
						TP_all ++;
					}
					else if((jv.getExpertClassification().equals(ExpertClassification.B) || jv.getExpertClassification().equals(ExpertClassification.LB)) && jv.getJudgment().getClassification().equals(Classification.Pathogn))
					{
						FP ++;
						FP_all ++;
					}
					else if((jv.getExpertClassification().equals(ExpertClassification.P) || jv.getExpertClassification().equals(ExpertClassification.LP)) && jv.getJudgment().getClassification().equals(Classification.Benign))
					{
						FN ++;
						FN_all ++;
					}
				}
			}
			
			System.out.println(mvl + "\t" + TN + "\t" + TP + "\t" + FP + "\t" + FN + "\t" + getTPR(TP, FN) + "\t" + getTNR(TN, FP)+ "\t" + getPPV(TP, FP)+ "\t" + getNPV(TN, FN));
			
			
		}
		
		System.out.println("TOTAL" + "\t" + TN_all + "\t" + TP_all + "\t" + FP_all + "\t" + FN_all + "\t" + getTPR(TP_all, FN_all) + "\t" + getTNR(TN_all, FP_all)+ "\t" + getPPV(TP_all, FP_all)+ "\t" + getNPV(TN_all, FN_all));

	}
	
	private static String getTPR(int TP, int FN)
	{
		return TP+FN == 0 ? "" : (int)Math.round((double)TP/(TP+FN)*100) + "%";
	}
	
	private static String getTNR(int TN, int FP)
	{
		return TN+FP == 0 ? "" : (int)Math.round((double)TN/(FP+TN)*100) + "%";
	}
	
	private static String getPPV(int TP, int FP)
	{
		return TP+FP == 0 ? "" : (int)Math.round((double)TP/(TP+FP)*100) + "%";
	}
	
	private static String getNPV(int TN, int FN)
	{
		return TN+FN == 0 ? "" : (int)Math.round((double)TN/(TN+FN)*100) + "%";
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
