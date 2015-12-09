package org.molgenis.calibratecadd.support;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.HashMap;
import java.util.Scanner;

import org.molgenis.data.annotation.joeri282exomes.Judgment;
import org.molgenis.data.annotation.joeri282exomes.Judgment.Classification;
import org.molgenis.data.annotation.joeri282exomes.Judgment.Method;

public class PONP2Results
{
	
	enum PonpClassification {
		Unknown, Neutral, Pathogenic
	}

	//chr:pos,ref,alt to PonpClassification
	HashMap<String, PonpClassification> ponpFile;
	
	/**
	 * 
	 * Example:
	 * 
	 * #Query	MappedGene	MappedVariation	ProbabilityOfPathogenicity	StandardError	Prediction	Annotations
	 * 12:32949167,T,C	ENSG00000057294	I789V	0.152	0.051	Neutral	 
	 * 2:220285283,C,G	ENSG00000175084	P268A	0.497	0.105	Unknown	 
	 * 18:28660261,C,T	ENSG00000134755	A441T	0.325	0.130	Unknown	 
	 * @throws FileNotFoundException 
	 * 
	 */
	public PONP2Results(File ponp2results) throws FileNotFoundException
	{
		ponpFile = new HashMap<String, PonpClassification>();
		Scanner s = new Scanner(ponp2results);
		s.nextLine();
		while(s.hasNextLine())
		{
			String line = s.nextLine();
			String[] linesplit = line.split("\t", -1);
			ponpFile.put(linesplit[0], PonpClassification.valueOf(linesplit[5]));
		}
		s.close();
	}
	
	public Judgment classifyVariantUsingPONP2Results(String chr, String pos, String ref, String alt) throws Exception
	{
		String key = chr+":"+pos+","+ref+","+alt;
		if(ponpFile.containsKey(key))
		{
			PonpClassification pc = ponpFile.get(key);
			if(pc.equals(PonpClassification.Unknown))
			{
				throw new VariantClassificationException("PONP2 result is 'unknown'a");
			}
			else if(pc.equals(PonpClassification.Neutral))
			{
				return new Judgment(Classification.Benign, Method.calibrated, "PONP2 result 'Neutral'");
			}
			else if(pc.equals(PonpClassification.Pathogenic))
			{
				return new Judgment(Classification.Pathogn, Method.calibrated, "PONP2 result 'Pathogenic'");
			}
			else
			{
				throw new VariantClassificationException("Unknown PONP2 result: " + pc);
			}
		}
		else
		{
			throw new VariantClassificationException("No PONP2 result");
		}
	}
	
	
}
