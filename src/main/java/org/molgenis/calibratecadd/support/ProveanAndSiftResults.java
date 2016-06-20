package org.molgenis.calibratecadd.support;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.HashMap;
import java.util.Scanner;

import org.molgenis.data.annotation.entity.impl.gavin.Judgment;
import org.molgenis.data.annotation.entity.impl.gavin.Judgment.Classification;
import org.molgenis.data.annotation.entity.impl.gavin.Judgment.Method;

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
			if(!linesplit[1].isEmpty())
			{
				ProveanClassification proveanPred = ProveanClassification.valueOf(linesplit[1]);
				proveanFile.put(linesplit[0], proveanPred);
			}
			if(!linesplit[2].isEmpty())
			{
				SiftClassification siftPred = SiftClassification.valueOf(linesplit[2]);
				siftFile.put(linesplit[0], siftPred);
			}
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
				return new Judgment(Classification.Benign, Method.calibrated, key, "PROVEAN result 'Neutral'");
			}
			else if(pc.equals(ProveanClassification.Deleterious))
			{
				return new Judgment(Classification.Pathogenic, Method.calibrated, key, "PROVEAN result 'Deleterious'");
			}
			else
			{
				return new Judgment(Classification.VOUS, Method.calibrated, key, "Unknown PROVEAN result: " + pc);
			}
		}
		else
		{
			return new Judgment(Classification.VOUS, Method.calibrated, key, "No PROVEAN result");
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
				return new Judgment(Classification.Benign, Method.calibrated, key, "SIFT result 'Tolerated'");
			}
			else if(pc.equals(SiftClassification.Damaging))
			{
				return new Judgment(Classification.Pathogenic, Method.calibrated, key, "SIFT result 'Damaging'");
			}
			else
			{
				return new Judgment(Classification.VOUS, Method.calibrated, key, "SIFT result '"+pc+"'");
			}
		}
		else
		{
			return new Judgment(Classification.VOUS, Method.calibrated, key, "No SIFT result");
		}
	}
	
	
}
