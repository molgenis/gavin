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
import org.molgenis.data.annotation.joeri282exomes.Judgment.Classification;
import org.molgenis.data.annotation.joeri282exomes.Judgment.Confidence;

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
		if(impact == null)
		{
			throw new InsufficientDataException("Missing impact for gene " + gene + " so we don't judge");
		}
		
		if(!geneToEntry.containsKey(gene))
		{
			return naiveClassifyVariant(gene, MAF, impact, CADDscore);
		//	throw new NoDataForGeneException("Cannot classify variant, no calibration data for " + gene);
		}
		
		if(geneToEntry.get(gene).PathoMAFThreshold == null)
		{
			throw new InsufficientDataException("Missing MAF for gene " + gene + " so we don't judge");
		}

		CCGGEntry entry = geneToEntry.get(gene);
		CCGGEntry.Category category = entry.category;
		
		/**
		 * High confidence tranche, strict settings
		 */
		
		// MAF filter
		if(MAF > 0 && MAF > (entry.PathoMAFThreshold*100))
		{
			return new Judgment(Classification.Benign, Confidence.High, "MAF of " + MAF + " greather than 100 times the pathogenic MAF "+(entry.PathoMAFThreshold+" "));
		}
		
		// "I-class genes"
		if(category.equals(Category.I1) && impact.equals(Impact.HIGH))
		{
			return new Judgment(Judgment.Classification.Pathogn,  Confidence.High, "Variant of HIGH impact, while no high impact variants known in population" + " (MAF " + MAF + ")");
		}
		else if(category.equals(Category.I2) && (impact.equals(Impact.MODERATE) || impact.equals(Impact.HIGH)))
		{
			return new Judgment(Judgment.Classification.Pathogn,  Confidence.High, "Variant of HIGH/MODERATE impact, while no high/moderate impact variants known in population" + " (MAF " + MAF + ")");
		}
		else if(category.equals(Category.I3) && (impact.equals(Impact.LOW) || impact.equals(Impact.MODERATE) || impact.equals(Impact.HIGH)))
		{
			return new Judgment(Judgment.Classification.Pathogn,  Confidence.High, "Variant of HIGH/MODERATE/LOW impact, while no high/moderate/low impact variants known in population" + " (MAF " + MAF + ")");
		}
		else if(impact.equals(Impact.MODIFIER))
		{
			return new Judgment(Judgment.Classification.Benign,  Confidence.High, "Variant impact of MODIFIER, unlikely to be pathogenic" + " (MAF " + MAF + ")");
		}

		// "C-class genes"
		if((category.equals(Category.C1) || category.equals(Category.C2)))
		{
			if(CADDscore != null && CADDscore > entry.MeanPathogenicCADDScore)
			{
				return new Judgment(Judgment.Classification.Pathogn,  Confidence.High, "CADD score of " + CADDscore + " higher than the mean pathogenic score of " + entry.MeanPathogenicCADDScore + " (and MAF " + MAF + ")");
			}
		}
		else if((category.equals(Category.C3) || category.equals(Category.C4) || category.equals(Category.C5)))
		{
			if(CADDscore != null && CADDscore > entry.Spec95thPerCADDThreshold)
			{
				return new Judgment(Judgment.Classification.Pathogn,  Confidence.High, "CADD score of " + CADDscore + " higher than the 95% specificity threhold " + entry.Spec95thPerCADDThreshold + " (and MAF " + MAF + ")");
			}
		}

		/**
		 * Medium confidence tranche, less strict settings
		 */
		
		if(MAF > 0 && MAF > (entry.PathoMAFThreshold))
		{
			return new Judgment(Classification.Benign, Confidence.Medium, "MAF > (entry.PathoMAFThreshold)");
		}
		
		if((category.equals(Category.C1) || category.equals(Category.C2)))
		{
			if(CADDscore != null && CADDscore > entry.MeanPopulationCADDScore)
			{
				return new Judgment(Judgment.Classification.Pathogn, Confidence.Medium, "CADDscore > entry.MeanPopulationCADDScore");
			}
		}
		else if((category.equals(Category.C3) || category.equals(Category.C4) || category.equals(Category.C5)))
		{
			if(CADDscore != null && CADDscore > entry.MeanPathogenicCADDScore)
			{
				return new Judgment(Judgment.Classification.Pathogn, Confidence.Medium, "CADDscore > entry.MeanPathogenicCADDScore");
			}
		}
		
		/**
		 * Low confidence tranche, least strict settings
		 */
		
		if(MAF > 0 && MAF > (entry.PathoMAFThreshold/100))
		{
			return new Judgment(Classification.Benign, Confidence.Low, "MAF > (entry.PathoMAFThreshold/100)");
		}
		
		if((category.equals(Category.C1) || category.equals(Category.C2)))
		{
			if(CADDscore != null && CADDscore > entry.Sens95thPerCADDThreshold)
			{
				return new Judgment(Judgment.Classification.Pathogn, Confidence.Low, "CADDscore > entry.Sens95thPerCADDThreshold");
			}
		}
		else if((category.equals(Category.C3) || category.equals(Category.C4) || category.equals(Category.C5)))
		{
			if(CADDscore != null && CADDscore > entry.MeanPopulationCADDScore)
			{
				return new Judgment(Judgment.Classification.Pathogn, Confidence.Low, "CADDscore > entry.MeanPopulationCADDScore");
			}
		}
	
	//	System.out.println("!! cant judge " + gene + " " + MAF + " " + impact + " " + CADDscore);
		
		return naiveClassifyVariant(gene, MAF, impact, CADDscore);
		
	//	throw new InsufficientDataException("No further calibration data for gene " + gene + " so we don't judge");
	}
	
	
	public Judgment naiveClassifyVariant(String gene, Double MAF, Impact impact, Double CADDscore) throws Exception
	{
		if(MAF > 0.00474)
		{
			return new Judgment(Judgment.Classification.Benign, Confidence.Naive, "MAF > 0.00474");
		}
//		if(impact.equals(impact.equals(Impact.MODIFIER) || impact.equals(Impact.LOW)))
//		{
//			return new Judgment(Judgment.Classification.Benign, Confidence.Naive, "Impact.MODIFIER || Impact.LOW");
//		}
		else
		{
			if(CADDscore != null && CADDscore > 34)
			{
				return new Judgment(Judgment.Classification.Pathogn, Confidence.Naive, "CADDscore > 34");
			}
			else if(CADDscore != null && CADDscore < 2)
			{
				return new Judgment(Judgment.Classification.Benign, Confidence.Naive, "CADDscore < 2");
			}
			else
			{
				throw new InsufficientDataException("Unable to naively classify " + gene + "! ");
			}
		}
	}
	
	public static Set<String> getGenesFromAnn(String ann) throws Exception
	{
		Set<String> genes = new HashSet<String>();
		String[] annSplit = ann.split(",", -1);
		for(String oneAnn : annSplit)
		{
			String[] fields = oneAnn.split("\\|", -1);
			String gene = fields[3];
			genes.add(gene);
		}
		if(genes.size() == 0)
		{
			throw new Exception("No genes for " + ann);
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
		String findAnn = getAnn(ann, gene, allele);
		if(findAnn == null)
		{
			System.out.println("FAILED TO GET IMPACT FOR " + ann);
			return null;
		}
		else
		{
			String[] fields = findAnn.split("\\|", -1);
			String impact = fields[2];
			return Impact.valueOf(impact);
		}
	}
	
	public static String getAnn(String ann, String gene, String allele) throws Exception
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
			return oneAnn;
		}
		System.out.println("warning: annotation could not be found for " + gene + ", allele=" + allele + ", ann=" + ann);
		return null;
	}

	public static void main(String[] args) throws Exception
	{
		File ccgg = new File(args[0]);
		new CCGGUtils(ccgg);

	}

}
