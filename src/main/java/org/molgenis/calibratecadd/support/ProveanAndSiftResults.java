package org.molgenis.calibratecadd.support;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.HashMap;
import java.util.Scanner;

import org.molgenis.data.annotation.joeri282exomes.Judgment;
import org.molgenis.data.annotation.joeri282exomes.Judgment.Classification;
import org.molgenis.data.annotation.joeri282exomes.Judgment.Method;

public class ProveanAndSiftResults
{
	
	enum ProveanClassification {
		Neutral, Deleterious
	}
	
	enum SiftClassification {
		Tolerated, Damaging, NA
	}

	//chr,pos,ref,alt to ProveanClassification
	HashMap<String, ProveanClassification> proveanFile;
	
	//chr,pos,ref,alt to SiftClassification
	HashMap<String, SiftClassification> siftFile;
	
	/**
	 * 
	 * Example:
	 * 
	 * INPUT	PROVEAN PREDICTION (cutoff=-2.5)	SIFT PREDICTION (cutoff=0.05)
	 * 1,154574443,C,T	Neutral	Tolerated
	 * 1,237947201,G,A	Neutral	Tolerated
	 * 2,220285283,C,G	Deleterious	Damaging
	 * 18,28660261,C,T	Deleterious	Damaging
	 * 
	 */
	public ProveanAndSiftResults(File proveanAndSiftResults) throws Exception
	{
		proveanFile = new HashMap<String, ProveanClassification>();
		siftFile = new HashMap<String, SiftClassification>();
		Scanner s = new Scanner(proveanAndSiftResults);
		s.nextLine(); //skip header
		while(s.hasNextLine())
		{
			String line = s.nextLine();
			String[] linesplit = line.split("\t", -1);
			ProveanClassification proveanPred = ProveanClassification.valueOf(linesplit[1]);
			SiftClassification siftPred = SiftClassification.valueOf(linesplit[2]);
			proveanFile.put(linesplit[0], proveanPred);
			siftFile.put(linesplit[0], siftPred);
		}
		s.close();
	}
	
	public Judgment classifyVariantUsingProveanResults(String chr, String pos, String ref, String alt) throws Exception
	{
		String key = chr+","+pos+","+ref+","+alt;
		if(proveanFile.containsKey(key))
		{
			ProveanClassification pc = proveanFile.get(key);
			if(pc.equals(ProveanClassification.Neutral))
			{
				return new Judgment(Classification.Benign, Method.calibrated, "PROVEAN result 'Neutral'");
			}
			else if(pc.equals(ProveanClassification.Deleterious))
			{
				return new Judgment(Classification.Pathogn, Method.calibrated, "PROVEAN result 'Deleterious'");
			}
			else
			{
				throw new VariantClassificationException("Unknown PROVEAN result: " + pc);
			}
		}
		else
		{
			throw new VariantClassificationException("No PROVEAN result");
		}
	}
	
	public Judgment classifyVariantUsingSiftResults(String chr, String pos, String ref, String alt) throws Exception
	{
		String key = chr+","+pos+","+ref+","+alt;
		if(siftFile.containsKey(key))
		{
			SiftClassification pc = siftFile.get(key);
			if(pc.equals(SiftClassification.Tolerated))
			{
				return new Judgment(Classification.Benign, Method.calibrated, "SIFT result 'Tolerated'");
			}
			else if(pc.equals(SiftClassification.Damaging))
			{
				return new Judgment(Classification.Pathogn, Method.calibrated, "SIFT result 'Damaging'");
			}
			else
			{
				throw new VariantClassificationException("Unknown SIFT result: " + pc);
			}
		}
		else
		{
			throw new VariantClassificationException("No SIFT result");
		}
	}
	
	
}
