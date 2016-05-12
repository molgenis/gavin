package org.molgenis.calibratecadd.support;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.HashMap;
import java.util.Scanner;

import org.molgenis.data.annotation.joeri282exomes.Judgment;
import org.molgenis.data.annotation.joeri282exomes.Judgment.Classification;
import org.molgenis.data.annotation.joeri282exomes.Judgment.Method;

public class MSCResults
{
	


	//genename to MSC threshold
	HashMap<String, Double> mscFile;
	
	/**
	 * 
	 */
	public MSCResults(File mscResults) throws FileNotFoundException
	{
		mscFile = new HashMap<String, Double>();
		Scanner s = new Scanner(mscResults);
		s.nextLine(); //skip "Gene	CADD_MSC"
		while(s.hasNextLine())
		{
			String line = s.nextLine();
			String[] linesplit = line.split("\t", -1);
			if(!linesplit[1].equals("NA"))
			{
				mscFile.put(linesplit[0], Double.parseDouble(linesplit[1]));
			}
			else
			{
				mscFile.put(linesplit[0], -1.0);
			}
		}
		s.close();
	}
	
	public Judgment classifyVariantUsingMSCResults(String geneName, Double caddScore) throws Exception
	{
		geneName = geneName.toUpperCase();
		if(caddScore == null)
		{
			return new Judgment(Classification.VOUS, Method.calibrated, geneName, "CADD score NULL for" + geneName);
		}
		else if(mscFile.containsKey(geneName) && mscFile.get(geneName) > 0)
		{
			if(caddScore > mscFile.get(geneName))
			{
				return new Judgment(Classification.Pathogn, Method.calibrated, geneName, "CADD score "+caddScore+" above MSC threshold " + mscFile.get(geneName));
			}
			else
			{
				return new Judgment(Classification.Benign, Method.calibrated, geneName, "CADD score "+caddScore+" below MSC threshold " + mscFile.get(geneName));
			}
		}
		else if(mscFile.containsKey(geneName) && mscFile.get(geneName) < 0)
		{
			return new Judgment(Classification.VOUS, Method.calibrated, geneName, "CADD score "+caddScore+", gene in MSC list but threshold unknown for " + geneName);
		}
		else
		{
			System.out.println("WARNING: no MSC threshold for gene " + geneName + ", try to add it?");
			return new Judgment(Classification.VOUS, Method.calibrated, geneName, "CADD score "+caddScore+" but no MSC threshold for gene " + geneName);
		}
	}
	
	
}
