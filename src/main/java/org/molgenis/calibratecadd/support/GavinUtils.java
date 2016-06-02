package org.molgenis.calibratecadd.support;

import java.io.File;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Scanner;
import java.util.Set;

import org.molgenis.data.Entity;

import org.molgenis.data.annotation.entity.impl.gavin.Judgment;
import org.molgenis.data.annotation.entity.impl.gavin.Judgment.Classification;
import org.molgenis.data.annotation.entity.impl.gavin.Judgment.Method;
import org.molgenis.data.annotation.entity.impl.snpEff.SnpEffRunner.Impact;
import org.molgenis.calibratecadd.support.CCGGEntry.Category;

public class GavinUtils
{
	HashMap<String, CCGGEntry> geneToEntry = new HashMap<String, CCGGEntry>();
	
	public GavinUtils(File ccgg) throws Exception
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

	public HashMap<String, CCGGEntry> getGeneToEntry()
	{
		return geneToEntry;
	}

	public Category getCategory(String gene)
	{
		return geneToEntry.get(gene).category;
	}
	
	public boolean contains(String gene)
	{
		return geneToEntry.containsKey(gene) ? true : false;
	}
	
	public Judgment classifyVariant(String gene, Double MAF, Impact impact, Double CADDscore)
	{
		//if we have no data for this gene, immediately fall back to the naive method
		if(!geneToEntry.containsKey(gene))
		{
			return naiveClassifyVariant(gene, MAF, impact, CADDscore);
		}

		CCGGEntry entry = geneToEntry.get(gene);
		CCGGEntry.Category category = entry.category;
		
		// MAF based classification, calibrated but slightly relaxed with a factor 10
		if(entry.PathoMAFThreshold != null && MAF > entry.PathoMAFThreshold * 10)
		{
			return new Judgment(Classification.Benign, Method.calibrated, gene, "Variant MAF of " + MAF + " is greater than the "+gene+" pathogenic 95th percentile MAF of "+ entry.PathoMAFThreshold + ".");
		}
		
		String mafReason = "the variant MAF of " + MAF + " is less than the "+gene+" pathogenic 95th percentile MAF of "+ entry.PathoMAFThreshold + ".";
		
		// Impact based classification, calibrated
		if(impact != null)
		{
			if(category.equals(Category.I1) && impact.equals(Impact.HIGH))
			{
				return new Judgment(Judgment.Classification.Pathogn,  Method.calibrated, gene, "Variant is of high impact, while there are no known high impact variants in the population. Also, " + mafReason);
			}
			else if(category.equals(Category.I2) && (impact.equals(Impact.MODERATE) || impact.equals(Impact.HIGH)))
			{
				return new Judgment(Judgment.Classification.Pathogn,  Method.calibrated, gene, "Variant is of high/moderate impact, while there are no known high/moderate impact variants in the population. Also, " + mafReason);
			}
			else if(category.equals(Category.I3) && (impact.equals(Impact.LOW) || impact.equals(Impact.MODERATE) || impact.equals(Impact.HIGH)))
			{
				return new Judgment(Judgment.Classification.Pathogn,  Method.calibrated, gene, "Variant is of high/moderate/low impact, while there are no known high/moderate/low impact variants in the population. Also, " + mafReason);
			}
			else if(impact.equals(Impact.MODIFIER))
			{
				return new Judgment(Judgment.Classification.Benign,  Method.calibrated, gene, "Variant is of 'modifier' impact, and therefore unlikely to be pathogenic. However, " + mafReason);
			}
		}

		// CADD score based classification, calibrated
		if(CADDscore != null)
		{
			if((category.equals(Category.C1) || category.equals(Category.C2)))
			{
				if(CADDscore > entry.MeanPathogenicCADDScore)
				{
					return new Judgment(Judgment.Classification.Pathogn,  Method.calibrated, gene, "Variant CADD score of " + CADDscore + " is greater than the mean pathogenic score of " + entry.MeanPathogenicCADDScore + " in a gene for which CADD scores are informative. Also, " + mafReason);
				}
//				else if(CADDscore < entry.MeanPopulationCADDScore)
//				{
//					return new Judgment(Judgment.Classification.Benign,  Method.calibrated, gene, "Variant CADD score of " + CADDscore + " is less than the mean population score of " + entry.MeanPathogenicCADDScore + " in a gene for which CADD scores are informative, although " + mafReason);
//				}
			}
			else if((category.equals(Category.C3) || category.equals(Category.C4) || category.equals(Category.C5)))
			{
				if(CADDscore > entry.Spec95thPerCADDThreshold)
				{
					return new Judgment(Judgment.Classification.Pathogn,  Method.calibrated, gene, "Variant CADD score of " + CADDscore + " is greater than the 95% specificity threhold of " + entry.Spec95thPerCADDThreshold + " for this gene. Also, " + mafReason);
				}
//				else if(CADDscore < entry.Sens95thPerCADDThreshold)
//				{
//					return new Judgment(Judgment.Classification.Benign,  Method.calibrated, gene, "Variant CADD score of " + CADDscore + " is less than the 95% sensitivity threhold of " + entry.MeanPathogenicCADDScore + " for this gene, although " + mafReason);
//				}
			}
		}
		
		//if everything so far has failed, we can still fall back to the naive method
		return naiveClassifyVariant(gene, MAF, impact, CADDscore);
	}
	
	
	public Judgment naiveClassifyVariant(String gene, Double MAF, Impact impact, Double CADDscore)
	{
		if(MAF > 0.00474)
		{
			return new Judgment(Judgment.Classification.Benign, Method.genomewide, gene, "Variant MAF of "+MAF+" is not rare enough to generally be considered pathogenic.");
		}
		if(impact.equals(Impact.MODIFIER))
		{
			return new Judgment(Judgment.Classification.Benign, Method.genomewide, gene, "Variant is of 'modifier' impact, and therefore unlikely to be pathogenic.");
		}
		else
		{
			if(CADDscore != null && CADDscore > 15)
			{
				return new Judgment(Judgment.Classification.Pathogn, Method.genomewide, gene, "Variant MAF of "+MAF+" is rare enough to be potentially pathogenic and the CADDscore of "+CADDscore+ " is greater than a global threshold of 15.");
			}
			else if(CADDscore != null && CADDscore <= 15)
			{
				return new Judgment(Judgment.Classification.Benign, Method.genomewide, gene, "CADDscore of "+CADDscore+ " is less than a global threshold of 15, although the variant MAF of "+MAF+" is rare enough to be potentially pathogenic.");
			}
			else
			{
				return new Judgment(Judgment.Classification.VOUS, Method.genomewide, gene, "Unable to classify variant as benign or pathogenic. The combination of "+impact+" impact, a CADD score "+CADDscore +" and MAF of " + MAF + " in " + gene + " is inconclusive.");
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
			System.out.println("WARNING: failed to get impact for gene '"+gene+"', allele '"+allele+"' in " + ann);
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
		System.out.println("WARNING: annotation could not be found for " + gene + ", allele=" + allele + ", ann=" + ann);
		return null;
	}

	public static void main(String[] args) throws Exception
	{
		File ccgg = new File(args[0]);
		new GavinUtils(ccgg);

	}

}
