package org.molgenis.data.annotation.joeri282exomes;

import java.io.File;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Scanner;
import java.util.Set;

import org.molgenis.calibratecadd.support.CaddScoreMissingException;
import org.molgenis.calibratecadd.support.InsufficientDataException;
import org.molgenis.calibratecadd.support.NoDataForGeneException;
import org.molgenis.data.Entity;
import org.molgenis.data.annotation.entity.impl.SnpEffAnnotator.Impact;
import org.molgenis.data.annotation.joeri282exomes.CCGGEntry.Category;

public class CCGGUtils
{
	HashMap<String, CCGGEntry> geneToEntry = new HashMap<String, CCGGEntry>();
	
	public CCGGUtils(File ccgg) throws Exception
	{	
		Scanner s = new Scanner(ccgg);
		
		//skip header
		s.nextLine();
		
		String line;
		while(s.hasNextLine())
		{
			line = s.nextLine();
			
			CCGGEntry e = new CCGGEntry(line);
			geneToEntry.put(e.gene, e);
		}
		
	}
	
	public Category getCategory(String gene)
	{
		return geneToEntry.get(gene).category;
	}
	
	public boolean contains(String gene)
	{
		return geneToEntry.containsKey(gene) ? true : false;
	}
	
	public Judgment classifyVariant(String gene, Double MAF, Impact impact, Double CADDscore) throws Exception
	{
		if(!geneToEntry.containsKey(gene))
		{
			throw new NoDataForGeneException("Cannot classify variant, no calibration data for " + gene);
		}
		
		CCGGEntry entry = geneToEntry.get(gene);
		CCGGEntry.Category category = entry.category;
		

		
		/**
		 * Apply MAF threshold first
		 * 100 is more strict in selecting (~1% FP), 10 is all right too (~4% FP)
		 */
		if(MAF > (entry.PathoMAFThreshold*100))
		{
			return new Judgment(Judgment.Classification.Likely_Benign, "todo");
		}
		

		/**
		 * "N-class genes"
		 */
		if((category.equals(Category.N1) || category.equals(Category.N2)))
		{
			
			
			throw new InsufficientDataException("Cannot classify variant, no futher calibration data for N1/N2 gene " + gene);
		}
		
		/**
		 * "T-class genes"
		 */
		if(category.equals(Category.T1) || category.equals(Category.T2))
		{			
			throw new InsufficientDataException("Cannot classify variant, no futher calibration data for T1/T2 gene " + gene);
		}
		
		
		
		
		/**
		 * "I-class genes"
		 */
		if(category.equals(Category.I1) && impact.equals(Impact.HIGH))
		{
			return new Judgment(Judgment.Classification.Pathogenic, "I1 and impact HIGH");
		}
		else if(category.equals(Category.I2) && (impact.equals(Impact.MODERATE) || impact.equals(Impact.HIGH)))
		{
			return new Judgment(Judgment.Classification.Pathogenic, "I2 and impact HIGH/MODERATE");
		}
		else if(category.equals(Category.I3) && (impact.equals(Impact.LOW) || impact.equals(Impact.MODERATE) || impact.equals(Impact.HIGH)))
		{
			return new Judgment(Judgment.Classification.Pathogenic, "I3 and impact HIGH/MODERATE/LOW");
		}
		else if(category.equals(Category.I1) || category.equals(Category.I2) || category.equals(Category.I3))
		{
			throw new InsufficientDataException("Cannot classify variant, no futher calibration data for I1/I2/I3 gene " + gene);
		}
		
		
		
		
		/**
		 * "C-class genes"
		 */
		if((category.equals(Category.C1) || category.equals(Category.C2)))
		{
			if(CADDscore == null)
			{
				throw new CaddScoreMissingException("Cannot classify variant, need to have CADD score for " + gene);
			}
			if(CADDscore > entry.Spec95thPerCADDThreshold)
			{
				return new Judgment(Judgment.Classification.Pathogenic, "CADD score " + CADDscore + " higher than 95% specificity threhold " + entry.Spec95thPerCADDThreshold);
			}
			else if(CADDscore > entry.MeanPathogenicCADDScore)
			{
				return new Judgment(Judgment.Classification.Likely_Pathogenic, "todo");
			}
			else if(CADDscore < entry.Sens95thPerCADDThreshold)
			{
				return new Judgment(Judgment.Classification.Benign, "todo");
			}
			else if(CADDscore < entry.MeanPopulationCADDScore)
			{
				return new Judgment(Judgment.Classification.Likely_Benign, "todo");
			}
			else
			{
				//justification?
			//	System.out.println("impact: " + impact);
				return new Judgment(Judgment.Classification.Likely_Pathogenic, "todo");
			}
		}
		
		if((category.equals(Category.C3) || category.equals(Category.C4)))
		{
			if(CADDscore == null)
			{
				throw new CaddScoreMissingException("Cannot classify variant, need to have CADD score for " + gene);
			}
			if(CADDscore > entry.Spec95thPerCADDThreshold)
			{
				return new Judgment(Judgment.Classification.Pathogenic, "CADD score " + CADDscore + " higher than 95% specificity threhold " + entry.Spec95thPerCADDThreshold);
			}
			else if(CADDscore > entry.MeanPathogenicCADDScore)
			{
				return new Judgment(Judgment.Classification.Likely_Pathogenic, "todo");
			}
			else if(CADDscore < entry.Sens95thPerCADDThreshold)
			{
				return new Judgment(Judgment.Classification.Benign, "todo");
			}
//			else if(CADDscore < entry.MeanPopulationCADDScore)
//			{
//				return new Judgment(Judgment.Classification.Likely_Benign, "todo");
//			}
		}
		
	//	if(1 == 1) { throw new NoDataForGeneException("QUITTING"); }
	
		
		throw new InsufficientDataException("Cannot classify variant in gene " + gene);
	}
	
	public static Set<String> getGenesFromAnn(String ann)
	{
		Set<String> genes = new HashSet<String>();
		String[] annSplit = ann.split(",", -1);
		for(String oneAnn : annSplit)
		{
			String[] fields = oneAnn.split("\\|", -1);
			String gene = fields[3];
			genes.add(gene);
		}
		return genes;
	}
	
	public static Double getInfoForAllele(Entity record, String infoField, String altAllele) throws Exception
	{
		String info_STR = record.get(infoField) == null ? null : record.get(infoField).toString();
		if(info_STR == null)
		{
			return null;
		}
		String[] alts = record.getString("ALT").split(",", -1);
		String[] info_split = info_STR.split(",", -1);
	
		if(alts.length != info_split.length)
		{
			throw new Exception("length of alts not equal to length of info field for " + record);
		}
		
		for (int i = 0; i < alts.length; i++)
		{
			if(alts[i].equals(altAllele))
			{
				return  (info_split[i] != null && !info_split[i].equals(".")) ? Double.parseDouble(info_split[i]) : null;
			}
		}
		return null;
	}
	
	public static Impact getImpact(String ann, String gene, String allele) throws Exception
	{
		String[] annSplit = ann.split(",", -1);
		for(String oneAnn : annSplit)
		{
			String[] fields = oneAnn.split("\\|", -1);
			String geneFromAnn = fields[3];
			if(!gene.equals(geneFromAnn))
			{
				continue;
			}
			String alleleFromAnn = fields[0];
			if(!allele.equals(alleleFromAnn))
			{
				continue;
			}
			String impact = fields[2];
			return Impact.valueOf(impact);
		}
		System.out.println("warning: impact could not be determined for " + gene + ", allele=" + allele + ", ann=" + ann);
		return null;
	}

	public static void main(String[] args) throws Exception
	{
		File ccgg = new File(args[0]);
		new CCGGUtils(ccgg);

	}

}
