package org.molgenis.calibratecadd;

import java.io.File;
import java.util.HashMap;
import java.util.Set;

import org.molgenis.calibratecadd.support.CaddScoreMissingException;
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

			for(String gene : genes)
			{
				if(ccgg.contains(gene))
				{
					Impact impact = CCGGUtils.getImpact(ann, gene, alt);
					System.out.println("going to classify: gene=" + gene + ", impact=" + impact + ", MAF=" + MAF + ", CADD=" + CADDscore);
					try{
					Judgment j = ccgg.classifyVariant(gene, MAF, impact, CADDscore);
					System.out.println("we say: " + j.getClassification() + ", mvl says: " + classification);
					addToMVLResults(j.getClassification(), classification, mvl);
					}
					catch(CaddScoreMissingException e)
					{
						System.out.println("unfortunately, missing cadd score to classify this variant");
					}
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


}
