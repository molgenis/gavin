package org.molgenis.data.annotation.joeri282exomes;

import java.io.File;
import java.util.HashMap;
import java.util.Scanner;

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
	
	public Judgment classifyVariant(String gene, Double MAF, Impact impact, Double CADDscore) throws Exception
	{
		//System.out.println("going to classify variant with gene="+gene+", maf=" + MAF + ", impact=" + impact + ", cadd=" + CADDscore);
		
		if(!geneToEntry.containsKey(gene))
		{
			throw new Exception("Cannot classify variant, no calibration data for " + gene);
		}
		CCGGEntry entry = geneToEntry.get(gene);
		CCGGEntry.Category category = entry.category;
		
		if((category.equals(Category.N1) || category.equals(Category.N2)))
		{
			return new Judgment(Judgment.Classification.VOUS, "todo");
		}
		
		if(MAF > entry.PathoMAFThreshold)
		{
	//		System.out.println("MAF > entry.PathoMAFThreshold: " + MAF + " > " + entry.PathoMAFThreshold);
			return new Judgment(Judgment.Classification.VOUS, "todo");
		}
		
		
		
		
	//	System.out.println("category: " + category);
		
		if((category.equals(Category.C1) || category.equals(Category.C2)))
		{
			if(CADDscore == null)
			{
				//TODO print where missing !!
		//		System.out.println("missing cadd!");
				return new Judgment(Judgment.Classification.VOUS, "todo");
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
				return new Judgment(Judgment.Classification.VOUS, "todo");
			}
		}
		if((category.equals(Category.C3) || category.equals(Category.C4)))
		{
			if(CADDscore == null)
			{
				return new Judgment(Judgment.Classification.VOUS, "todo");
			}
			if(CADDscore > entry.Spec95thPerCADDThreshold)
			{
				//LP ?
				return new Judgment(Judgment.Classification.VOUS, "todo");
			}
			else if(CADDscore < entry.Sens95thPerCADDThreshold)
			{
				return new Judgment(Judgment.Classification.VOUS, "todo");
			}
			else
			{
				return new Judgment(Judgment.Classification.VOUS, "todo");
			}
		}
		else if(category.equals(Category.I1) && impact.equals(Impact.HIGH))
		{
			return new Judgment(Judgment.Classification.Pathogenic, "category.equals(Category.I1) && impact.equals(Impact.HIGH)) so Pathogenic");
		}
		else if(category.equals(Category.I2) && impact.equals(Impact.MODERATE))
		{
			return new Judgment(Judgment.Classification.Pathogenic, "category.equals(Category.I2) && impact.equals(Impact.MODERATE) so Pathogenic");
		}
		else if(category.equals(Category.I3) && impact.equals(Impact.LOW))
		{
			return new Judgment(Judgment.Classification.Pathogenic, "category.equals(Category.I3) && impact.equals(Impact.LOW) so Pathogenic");
		}
		
		return new Judgment(Judgment.Classification.VOUS, "todo");
	}
	
	

	public static void main(String[] args) throws Exception
	{
		File ccgg = new File(args[0]);
		new CCGGUtils(ccgg);

	}

}
