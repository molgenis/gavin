package org.molgenis.calibratecadd.support;

import java.nio.file.Files;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.util.HashMap;
import java.util.List;

import org.molgenis.calibratecadd.support.JudgedVariant.ExpertClassification;
import org.molgenis.data.annotation.joeri282exomes.Judgment.Classification;
import org.molgenis.data.annotation.joeri282exomes.Judgment.Method;

public class ProcessJudgedVariantMVLResults
{
	
	private static Integer grandTotalExpertClassified;
	private static Integer grandTotalOursClassified;
	private static Double MCCofAllMVLs;
	private static Integer _TN_all, _TP_all, _FP_all, _FN_all, _vousB, _vousP;
	
	public static void printResults(HashMap<String, List<JudgedVariant>> judgedMVLVariants, String toolName, String datasetName) throws Exception
	{
		printCountsOfExpertMVLClassifications(judgedMVLVariants);
		printCountsOfCCGGMVLClassifications(judgedMVLVariants, null);
		calculateAndPrint_FP_FN_stats(judgedMVLVariants, null);
		printYield();
		printCountsOfCCGGMVLClassifications(judgedMVLVariants, Method.calibrated);
		calculateAndPrint_FP_FN_stats(judgedMVLVariants, Method.calibrated);
		printCountsOfCCGGMVLClassifications(judgedMVLVariants, Method.naive);
		calculateAndPrint_FP_FN_stats(judgedMVLVariants, Method.naive);
		reportVOUScounts(judgedMVLVariants, Method.calibrated);
		printVOUSresults(judgedMVLVariants, Method.calibrated);
		reportVOUScounts(judgedMVLVariants, Method.naive);
		printVOUSresults(judgedMVLVariants, Method.naive);
		printFalseResults(judgedMVLVariants, Method.calibrated);
		printFalseResults(judgedMVLVariants, Method.naive);
		//printCodeForPieChart(toolName, datasetName);
		printCodeForDF(toolName, datasetName);
	}
	
	public static void printCodeForDF(String toolName, String datasetName) throws Exception
	{
		if(_TN_all == null || _TP_all == null || _FP_all == null || _FN_all == null || _vousB == null || _vousP == null || MCCofAllMVLs == null)
		{
			throw new Exception("Can only print R code when we have TN_all, TP_all, FP_all, FN_all, grandTotalExpertClassified, MCCofAllMVLs");
		}
		
		String toR = "row <- data.frame(Tool = \""+toolName+"\", Data = \""+datasetName+"\", MCC = "+MCCofAllMVLs+", TN = "+_TN_all+", TP = "+_TP_all+", FP = "+_FP_all+", FN = "+_FN_all+", VB = "+_vousB+", VP = "+_vousP+"); df <- rbind(df, row)\n";
		Files.write(Paths.get("/Users/jvelde/dfrows.r"), toR.getBytes(), StandardOpenOption.APPEND);
	}
	public static void printCodeForPieChart(String toolName, String datasetName) throws Exception
	{
		if(_TN_all == null || _TP_all == null || _FP_all == null || _FN_all == null || grandTotalExpertClassified == null)
		{
			throw new Exception("Can only print piechart R code when we have TN_all, TP_all, FP_all, FN_all, grandTotalExpertClassified");
		}
		
		//System.out.println("TOTAL <- " + grandTotalExpertClassified);
		System.out.println("TN <- "+_TN_all+"; TP <- "+_TP_all+"; FP <- "+_FP_all+"; FN <- "+_FN_all+"");
		System.out.println("NC <- TOTAL-TN-TP-FP-FN");
		System.out.println("slices <- c(TN, TP, FP, FN, NC)");
		System.out.println("lbls <- c(\"TN\",\"TP\",\"FP\",\"FN\", \"Not cl.\")");
		System.out.println("pct <- round(slices/sum(slices)*100)");
		System.out.println("lbls <- paste(lbls, pct)");
		System.out.println("lbls <- paste(lbls,\"%\",sep=\"\")");
		System.out.println("png(filename=\"~/"+toolName+"_"+datasetName+".png\",res=100, width=600, height=600)");
		System.out.println("pie(slices,labels = lbls, col=c(\"#008800\",\"#00ff00\",\"#ff0000\",\"#880000\",\"#cccccc\"), main=\"Applied "+toolName+" on "+datasetName+" (total "+grandTotalExpertClassified+" variants)\")");
		System.out.println("dev.off() ;");
	}
	
	//TODO is this correct? grandTotalExpertClassified may contain VOUS variants, which we cannot classify 'wrong' or 'right'.... hmm....
	public static void printYield() throws Exception
	{
		if(grandTotalExpertClassified == null || grandTotalOursClassified == null || MCCofAllMVLs == null)
		{
			throw new Exception("Can only calculate yield when we have grandTotalExpertClassified, grandTotalOursClassified, MCCofAllMVLs");
		}
		else
		{
			System.out.println("\nYield:\n" + grandTotalOursClassified + " / " + grandTotalExpertClassified + " * " + MCCofAllMVLs + " = " + ((double)grandTotalOursClassified/grandTotalExpertClassified*MCCofAllMVLs));
		}
	}
	
	public static void printFalseResults(HashMap<String, List<JudgedVariant>> judgedMVLVariants, Method method)
	{
		System.out.println("\nFalse hits, method: " + method);
		
		for(String mvl : judgedMVLVariants.keySet())
		{
			StringBuffer FN = new StringBuffer();
			
			for(JudgedVariant jv : judgedMVLVariants.get(mvl))
			{
				if(jv.getJudgment() != null && jv.getJudgment().getClassification().equals(Classification.Benign) && (jv.getExpertClassification().equals(ExpertClassification.P) || jv.getExpertClassification().equals(ExpertClassification.LP)) && jv.getJudgment().getConfidence().equals(method))
				{
					FN.append(jv.printVariant() + "\n");
				}
			}
			if(FN.length() > 0)
			{
				System.out.println(mvl + ", false negatives:" + "\n" + FN.toString());
			}
			
			
			StringBuffer FP = new StringBuffer();
			
			for(JudgedVariant jv : judgedMVLVariants.get(mvl))
			{
				if(jv.getJudgment() != null && jv.getJudgment().getClassification().equals(Classification.Pathogn) && (jv.getExpertClassification().equals(ExpertClassification.B) || jv.getExpertClassification().equals(ExpertClassification.LB)) && jv.getJudgment().getConfidence().equals(method))
				{
					FP.append(jv.printVariant() + "\n");
				}
			}
			if(FP.length() > 0)
			{
				System.out.println(mvl + ", false positives:" + "\n" + FP.toString());
			}
		}
	}
	
	public static void printVOUSresults(HashMap<String, List<JudgedVariant>> judgedMVLVariants, Method method)
	{
		System.out.println("\nVOUS variants, method: " + method);
		
		for(String mvl : judgedMVLVariants.keySet())
		{
			StringBuffer benignVous = new StringBuffer();
			
			for(JudgedVariant jv : judgedMVLVariants.get(mvl))
			{
				if(jv.getJudgment() != null && jv.getJudgment().getClassification().equals(Classification.Benign) && jv.getExpertClassification().equals(ExpertClassification.V) && jv.getJudgment().getConfidence().equals(method))
				{
					benignVous.append(jv.printVariant() + "\n");
				}
			}
			if(benignVous.length() > 0)
			{
				System.out.println(mvl + ", benign:" + "\n" + benignVous.toString());
			}
			
			
			StringBuffer pathoVous = new StringBuffer();
			
			for(JudgedVariant jv : judgedMVLVariants.get(mvl))
			{
				if(jv.getJudgment() != null && jv.getJudgment().getClassification().equals(Classification.Pathogn) && jv.getExpertClassification().equals(ExpertClassification.V) && jv.getJudgment().getConfidence().equals(method))
				{
					pathoVous.append(jv.printVariant() + "\n");
				}
			}
			if(pathoVous.length() > 0)
			{
				System.out.println(mvl + ", pathogenic:" + "\n" + pathoVous.toString());
			}
		}
	}
	
	public static void reportVOUScounts(HashMap<String, List<JudgedVariant>> judgedMVLVariants, Method method)
	{
		System.out.println("\nClassifications of VOUS variants, method: " + method);
		
		System.out.println("\t" + "Benign" + "\t" + "Pathogn");
		
		int totalBenign = 0;
		int totalPathogn = 0;
		for(String mvl : judgedMVLVariants.keySet())
		{
			int benignVOUSforMVL = 0;
			int pathogenicVOUSforMVL = 0;
			for(JudgedVariant jv : judgedMVLVariants.get(mvl))
			{
				if(jv.getJudgment() != null && jv.getExpertClassification().equals(ExpertClassification.V) && jv.getJudgment().getConfidence().equals(method))
				{
					if(jv.getJudgment().getClassification().equals(Classification.Benign))
					{
						benignVOUSforMVL++;
						totalBenign++;
					}
					else if(jv.getJudgment().getClassification().equals(Classification.Pathogn))
					{
						pathogenicVOUSforMVL++;
						totalPathogn++;
					}
				}
			}
			System.out.println(mvl + "\t" + benignVOUSforMVL + "\t" + pathogenicVOUSforMVL);
		}
		System.out.println("TOTAL" + "\t" + totalBenign + "\t" + totalPathogn);
		
		
	}
	
	public static void calculateAndPrint_FP_FN_stats(HashMap<String, List<JudgedVariant>> judgedMVLVariants, Method method)
	{
		System.out.println("\nFalse posities & false negatives, method: " + (method == null ? "all" : method));
		int grandTotal = 0;
		
		System.out.println("\t" + "#TN" + "\t" + "#TP" + "\t" + "#FP" + "\t" + "#FN" + "\t" + "TPR" + "\t" + "TNR" + "\t" + "PPV" + "\t" + "NPV" + "\t" + "ACC" + "\t" + "MCC" + "\t" + "OPM");
		
		int TN_all = 0;
		int TP_all = 0;
		int FP_all = 0;
		int FN_all = 0;
		int vousB_all = 0;
		int vousP_all = 0;
		
		for(String mvl : judgedMVLVariants.keySet())
		{
			int TN = 0;
			int TP = 0;
			int FP = 0;
			int FN = 0;
			for(JudgedVariant jv : judgedMVLVariants.get(mvl))
			{
				if(jv.getJudgment() != null && (method == null || jv.getJudgment().getConfidence().equals(method)))
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
					else if((jv.getExpertClassification().equals(ExpertClassification.P) || jv.getExpertClassification().equals(ExpertClassification.LP)) && jv.getJudgment().getClassification().equals(Classification.VOUS))
					{
						vousP_all ++;
					}
					else if((jv.getExpertClassification().equals(ExpertClassification.B) || jv.getExpertClassification().equals(ExpertClassification.LB)) && jv.getJudgment().getClassification().equals(Classification.VOUS))
					{
						vousB_all ++;
					}
				}
			}
			
			System.out.println(mvl + "\t" + TN + "\t" + TP + "\t" + FP + "\t" + FN + "\t" + getTPR(TP, FN) + "\t" + getTNR(TN, FP)+ "\t" + getPPV(TP, FP) + "\t" + getNPV(TN, FN) + "\t" + getAcc(TP, TN, FP, FN) + "\t" + getMCC(TP, TN, FP, FN) + "\t" + getOPM(TP, TN, FP, FN));
			
			
		}
		
		System.out.println("TOTAL" + "\t" + TN_all + "\t" + TP_all + "\t" + FP_all + "\t" + FN_all + "\t" + getTPR(TP_all, FN_all) + "\t" + getTNR(TN_all, FP_all)+ "\t" + getPPV(TP_all, FP_all)+ "\t" + getNPV(TN_all, FN_all) + "\t" + getAcc(TP_all, TN_all, FP_all, FN_all) + "\t" + getMCC(TP_all, TN_all, FP_all, FN_all) + "\t" + getOPM(TP_all, TN_all, FP_all, FN_all));
		
		if(method == null) //mode = 'all'
		{
			MCCofAllMVLs = getMCC(TP_all, TN_all, FP_all, FN_all);
			_TN_all = TN_all;
			_TP_all = TP_all;
			_FP_all = FP_all;
			_FN_all = FN_all;
			_vousB = vousB_all;
			_vousP = vousP_all;
		}
		
	}
	
	private static String getTPR(int TP, int FN)
	{
		return TP+FN == 0 ? "-" : (int)Math.round((double)TP/(TP+FN)*100) + "%";
	}
	
	private static String getTNR(int TN, int FP)
	{
		return TN+FP == 0 ? "-" : (int)Math.round((double)TN/(FP+TN)*100) + "%";
	}
	
	private static String getPPV(int TP, int FP)
	{
		return TP+FP == 0 ? "-" : (int)Math.round((double)TP/(TP+FP)*100) + "%";
	}
	
	private static String getNPV(int TN, int FN)
	{
		return TN+FN == 0 ? "-" : (int)Math.round((double)TN/(TN+FN)*100) + "%";
	}
	
	private static String getAcc(int TP, int TN, int FP, int FN)
	{
		return TP+TN+FP+FN == 0 ? "-" : (int)Math.round((double)(TP+TN)/(TP+TN+FP+FN)*100) + "%";
	}
	
	// https://en.wikipedia.org/wiki/Matthews_correlation_coefficient
	private static double getMCC(int TP, int TN, int FP, int FN)
	{
		double aboveDiv = TP*TN-FP*FN;
		double belowDiv = Math.sqrt((double)(TP+FP)*(TP+FN)*(TN+FP)*(TN+FN));
		if(belowDiv == 0)
		{
			return 0; // arbitrarily set to 0, though could also return "-" or NULL or so..
		}
		double mcc = (aboveDiv/belowDiv);
		//double nMcc = ( 1 + mcc ) / 2;
		return mcc;
	}
	
	private static String MCCtoString(double MCC)
	{
		if(MCC == 0)
		{
			return "-";
		}
		else
		{
			return (int)Math.round((MCC)*100) + "%";
		}
	}
	
	private static String getOPM(int TP, int TN, int FP, int FN)
	{
		if(TP+FP == 0 || TP+FN == 0 || TN+FP == 0 || TN+FN == 0)
		{
			return "-";
		}
		double ppv =  (double)TP/(TP+FP);
		double npv = (double)TN/(TN+FN);
		double sens = (double)TP/(TP+FN); //same thing as TPR
		double spec = (double)TN/(FP+TN); //same thing as TNR
		double acc = (double)(TP+TN)/(TP+TN+FP+FN);
		
		double aboveDiv = TP*TN-FP*FN;
		double belowDiv = Math.sqrt((double)(TP+FP)*(TP+FN)*(TN+FP)*(TN+FN));
		double mcc = aboveDiv/belowDiv;
		double nMcc = ( 1 + mcc ) / 2;
//		System.out.println("("+ppv +"+"+ npv+") * ("+sens +"+"+ spec+") * ("+acc +"+"+ nMcc+")");
		double OPM = ( (ppv + npv) * (sens + spec) * (acc + nMcc) ) / 8;
		return (int)Math.round((OPM)*100) + "%";
	}

	public static void printCountsOfCCGGMVLClassifications(HashMap<String, List<JudgedVariant>> judgedMVLVariants, Method method)
	{
		System.out.println("\nCounts, method: " + (method == null ? "all" : method));
		int grandTotal = 0;
		int totalBenign = 0;
		int totalPathogenic = 0;
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
					if(jv.getJudgment().getClassification().equals(cl) && (method == null || jv.getJudgment().getConfidence().equals(method)))
					{
						if(jv.getJudgment().getClassification().equals(Classification.Benign))
						{
							totalBenign++;
						}
						if(jv.getJudgment().getClassification().equals(Classification.Pathogn))
						{
							totalPathogenic++;
						}
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
				if(method == null || jv.getJudgment().getConfidence().equals(method))
				{
					count++;
				}
			}
			System.out.print("\t" + count);
			count = 0;
		}
		System.out.println("\t" + grandTotal);
		grandTotalOursClassified = totalBenign + totalPathogenic;
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
		grandTotalExpertClassified = grandTotal;
	}

}
