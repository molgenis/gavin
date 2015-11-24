package org.molgenis.calibratecadd;

import java.io.File;
import java.util.HashMap;
import java.util.Set;

import org.molgenis.calibratecadd.support.CCGGException;
import org.molgenis.calibratecadd.support.MVLResults;
import org.molgenis.data.Entity;
import org.molgenis.data.annotation.entity.impl.SnpEffAnnotator.Impact;
import org.molgenis.data.annotation.joeri282exomes.CCGGUtils;
import org.molgenis.data.annotation.joeri282exomes.Judgment;
import org.molgenis.data.annotation.joeri282exomes.Judgment.Classification;
import org.molgenis.data.vcf.VcfRepository;

public class Step9_Validation
{
	public static void main(String[] args) throws Exception
	{
		new Step9_Validation(args[0], args[1]);
	}
	
	CCGGUtils ccgg;
	HashMap<String, MVLResults> mvlResults = new HashMap<String, MVLResults>();

	/**
	 * Take MVL (annotated with CADD, ExAC, SnpEff)
	 * Take CCGG thresholds
	 * Check if classification matches
	 * @throws Exception 
	 */
	public Step9_Validation(String ccggLoc, String mvlLoc) throws Exception
	{
		loadCCGG(ccggLoc);
		scanMVL(mvlLoc);
		reportResults();
		calculateOverallPerformance();
	}
	
	public void scanMVL(String mvlLoc) throws Exception
	{
		File mvlFile = new File(mvlLoc);
		if(!mvlFile.exists())
		{
			throw new Exception("MVL file "+mvlFile+" does not exist or is directory");
		}
		VcfRepository vcfRepo = new VcfRepository(mvlFile, "mvl");
		
		java.util.Iterator<Entity> vcfRepoIter = vcfRepo.iterator();

		while (vcfRepoIter.hasNext())
		{
			Entity record = vcfRepoIter.next();
			//System.out.println(record.toString());
			
			String alt = record.getString("ALT");
			if(alt.contains(","))
			{
				throw new Exception("Did not expect multiple alt alleles! " + record.toString());
			}
			
			Double getMAF = CCGGUtils.getInfoForAllele(record, "EXAC_AF", alt);
			double MAF = getMAF == null ? 0 : getMAF;
			
			Double CADDscore = CCGGUtils.getInfoForAllele(record, "CADD_SCALED", alt);
			
			String ann = record.getString("ANN");
			Set<String> genes = CCGGUtils.getGenesFromAnn(ann);
			
			String classification = record.getString("CLSF");
			String mvl = record.getString("MVL");
			
//			String id = record.getString("ID");
//			if( nc.newClsf.containsKey(id))
//			{
//				System.out.println((classification.equals(nc.newClsf.get(id))?"" : "---> !! DIFFERENCE: ") + id + " is found in VCF classified as " + classification + ", according to re-classified list: " + nc.newClsf.get(id));
//			}
			
			for(String gene : genes)
			{
		
		//		System.out.println("MVL: " + mvl + ", gene: " + gene + ", cat:" + ccgg.getCategory(gene));
				Impact impact = CCGGUtils.getImpact(ann, gene, alt);
		//		System.out.println("going to classify: gene=" + gene + ", impact=" + impact + ", MAF=" + MAF + ", CADD=" + CADDscore);
				
				try{
					Judgment j = ccgg.classifyVariant(gene, MAF, impact, CADDscore);
				//	System.out.println("we say: " + j.getClassification() + ", mvl says: " + classification);
					addToMVLResults(j.getClassification(), classification, mvl);
					}
				catch(CCGGException e)
				{
		//			System.out.println(e);
					addToMVLResults(null, classification, mvl);
				}
			}
		}
	}
	
	private void addToMVLResults(Classification ccggClasfc, String mvlClasfc, String mvl)
	{
		
		if(!mvlResults.containsKey(mvl))
		{
			MVLResults mvlR = new MVLResults();
			mvlResults.put(mvl, mvlR);
		}
		
		MVLResults mvlR = mvlResults.get(mvl);
		
		if(mvlClasfc.equals("B") || mvlClasfc.equals("LB"))
		{
			mvlR.nrOf_B_LB++;
			if(ccggClasfc == Classification.Benign || ccggClasfc == Classification.Likely_Benign)
			{
				mvlR.nrOf_B_LB_judgedAs_B_LB++;
			}
			if(ccggClasfc == Classification.Pathogenic || ccggClasfc == Classification.Likely_Pathogenic)
			{
				//"false positives"
				mvlR.nrOf_B_LB_judgedAs_P_LP++;
			}
		}
		
		if(mvlClasfc.equals("P") || mvlClasfc.equals("LP"))
		{
			mvlR.nrOf_P_LP++;
			if(ccggClasfc == Classification.Benign || ccggClasfc == Classification.Likely_Benign)
			{
				//"false negatives"
				mvlR.nrOf_P_LP_judgedAs_B_LB++;
			}
			if(ccggClasfc == Classification.Pathogenic || ccggClasfc == Classification.Likely_Pathogenic)
			{
				mvlR.nrOf_P_LP_judgedAs_P_LP++;
			}
		}
		

		if(mvlClasfc.equals("V"))
		{
			mvlR.nrOf_VOUS++;
			if(ccggClasfc == Classification.Benign || ccggClasfc == Classification.Likely_Benign)
			{
				mvlR.nrOf_VOUS_judgedAs_B_LB++;
			}
			if(ccggClasfc == Classification.Pathogenic || ccggClasfc == Classification.Likely_Pathogenic)
			{
				mvlR.nrOf_VOUS_judgedAs_P_LP++;
			}
		}
	}
	
	public void loadCCGG(String ccggLoc) throws Exception
	{
		File ccggFile = new File(ccggLoc);
		if(!ccggFile.exists())
		{
			throw new Exception("CCGG file "+ccggFile+" does not exist or is directory");
		}
		CCGGUtils ccggUtils = new CCGGUtils(ccggFile);
		this.ccgg = ccggUtils;
	}
	
	public void reportResults()
	{
		for(String mvl : mvlResults.keySet())
		{
			System.out.println("MVL: " + mvl);
			System.out.println(mvlResults.get(mvl).toString());
			System.out.println("FP: " +(mvlResults.get(mvl).getPercOfFalsePositives()));
			System.out.println("FN: " +(mvlResults.get(mvl).getPercOfFalseNegatives()));
			
		}
		
	}
	
	public void calculateOverallPerformance()
	{
		int correctlyClassifiedBenign = 0;
		int correctlyClassifiedPathogenic = 0;
		int wronglyClassifiedBenign = 0;
		int wronglyClassifiedPathogenic = 0;
		int nrOf_B_LB_P_LP_CCGGclassifiedVariants_acrossAllMVLs = 0;
		
		int nrOf_B_LB_CCGGclassifiedVariants_acrossAllMVLs = 0;
		int nrOf_P_LP_CCGGclassifiedVariants_acrossAllMVLs = 0;
		
		int nrOf_VOUS_CCGGclassifiedVariants_acrossAllMVLs = 0;
		int nrOf_VOUS_B_LB_CCGGclassifiedVariants_acrossAllMVLs = 0;
		int nrOf_VOUS_P_LP_CCGGclassifiedVariants_acrossAllMVLs = 0;
		
		int nrOfMVL_B_LB_P_LP_variants = 0;
		int nrOfMVL_VOUS_variants = 0;
		int nrOfMVL_P_LP_variants = 0;
		int nrOfMVL_B_LB_variants = 0;
		for(String mvl : mvlResults.keySet())
		{
			nrOfMVL_B_LB_P_LP_variants += mvlResults.get(mvl).nrOf_B_LB + mvlResults.get(mvl).nrOf_P_LP;
			nrOfMVL_P_LP_variants += mvlResults.get(mvl).nrOf_P_LP;
			nrOfMVL_B_LB_variants += mvlResults.get(mvl).nrOf_B_LB;
			
			nrOf_B_LB_P_LP_CCGGclassifiedVariants_acrossAllMVLs += mvlResults.get(mvl).nrOf_B_LB_P_LP_CCGGClassifiedVariants();
			nrOf_B_LB_CCGGclassifiedVariants_acrossAllMVLs += mvlResults.get(mvl).nrOf_B_LB_judgedAs_B_LB + mvlResults.get(mvl).nrOf_P_LP_judgedAs_B_LB;
			nrOf_P_LP_CCGGclassifiedVariants_acrossAllMVLs += mvlResults.get(mvl).nrOf_P_LP_judgedAs_P_LP + mvlResults.get(mvl).nrOf_B_LB_judgedAs_P_LP;
			
			nrOfMVL_VOUS_variants += mvlResults.get(mvl).nrOf_VOUS;		
			nrOf_VOUS_CCGGclassifiedVariants_acrossAllMVLs += mvlResults.get(mvl).nrOf_VOUS_CCGGClassifiedVariants();
			nrOf_VOUS_B_LB_CCGGclassifiedVariants_acrossAllMVLs += mvlResults.get(mvl).nrOf_VOUS_judgedAs_B_LB;
			nrOf_VOUS_P_LP_CCGGclassifiedVariants_acrossAllMVLs += mvlResults.get(mvl).nrOf_VOUS_judgedAs_P_LP;
			
			correctlyClassifiedBenign += mvlResults.get(mvl).nrOf_B_LB_judgedAs_B_LB;
			correctlyClassifiedPathogenic += mvlResults.get(mvl).nrOf_P_LP_judgedAs_P_LP;
			wronglyClassifiedBenign += mvlResults.get(mvl).nrOf_B_LB_judgedAs_P_LP;
			wronglyClassifiedPathogenic += mvlResults.get(mvl).nrOf_P_LP_judgedAs_B_LB;
			
		}

		
		System.out.println("\n\nTotal amount of B/LB/LP/P variants in original MVLs = " + nrOfMVL_B_LB_P_LP_variants + " ("+nrOfMVL_B_LB_variants+" benign, "+nrOfMVL_P_LP_variants+" pathogenic)");
		System.out.println("Total amount of B/LB/LP/P variants that CCGG-classified as B/LB/LP/P = " + nrOf_B_LB_P_LP_CCGGclassifiedVariants_acrossAllMVLs);
		
		System.out.println("\nCorrectly classified benign = " + correctlyClassifiedBenign);
		System.out.println("Wrongly classified benign (LB/B judged as LP/P, false positives) = " + wronglyClassifiedBenign);
		System.out.println("Correctly classified pathogenic = " + correctlyClassifiedPathogenic);
		System.out.println("Wrongly classified pathogenic (LP/P judged as LB/B, false negatives) = " + wronglyClassifiedPathogenic);

		System.out.println("\nAverage false positive rate in CCGG-classifications across all MVLs = " + Math.round(((double)wronglyClassifiedBenign/(double)(wronglyClassifiedBenign+correctlyClassifiedBenign))*100) + "%");
		System.out.println("Average false negative rate in CCGG-classifications across all MVLs = " + Math.round(((double)wronglyClassifiedPathogenic/(double)(wronglyClassifiedPathogenic+correctlyClassifiedPathogenic))*100) + "%");

		System.out.println("\nClassification of VOUS variants:");
		System.out.println("Total amount of VOUS variants in original MVLs = " + nrOfMVL_VOUS_variants);
		System.out.println("Total amount of VOUS variants that CCGG-classified as B/LB/LP/P = " + nrOf_VOUS_CCGGclassifiedVariants_acrossAllMVLs + " (" + nrOf_VOUS_B_LB_CCGGclassifiedVariants_acrossAllMVLs + " benign, " + nrOf_VOUS_P_LP_CCGGclassifiedVariants_acrossAllMVLs + " pathogenic)");

		
	}


}
