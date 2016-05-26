package org.molgenis.calibratecadd.support;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.HashMap;
import java.util.Scanner;

import org.molgenis.data.annotation.joeri282exomes.Judgment;
import org.molgenis.data.annotation.joeri282exomes.Judgment.Classification;
import org.molgenis.data.annotation.joeri282exomes.Judgment.Method;

public class MutationTaster2Results
{
	
	enum MutationTasterClassification {
		disease_causing, polymorphism
	}

	//chr:pos,ref,alt to MutationTasterClassification
	HashMap<String, MutationTasterClassification> muttasterFile;
	
	/**
	 * 
	 * Example:
	 * 
	 * chromosome	position	pred_index	allele_ref	allele_alt
	 * 18	29178610	disease_causing_automatic	C	T
	 * 18	29178554	disease_causing	C	T
	 * 10	92679011	disease_causing_automatic		T
	 * 
	 */
	public MutationTaster2Results(File muttaster2results) throws Exception
	{
		muttasterFile = new HashMap<String, MutationTasterClassification>();
		Scanner s = new Scanner(muttaster2results);
		s.nextLine(); //skip header
		while(s.hasNextLine())
		{
			String line = s.nextLine();
			String[] linesplit = line.split("\t", -1);
			
			//2	220285298	disease_causing_automatic	G	CG
			//2	220284876	disease_causing	C	T
			//2	220285088	polymorphism	C	T
			//2	220286064	polymorphism_automatic	C	T
			String chrPosRefAlt = linesplit[0] + ":" + linesplit[1] + "," + linesplit[3] + "," + linesplit[4];
			MutationTasterClassification prediction = null;
			if(linesplit[2].startsWith("disease_causing"))
			{
				prediction = MutationTasterClassification.disease_causing;
			}
			else if(linesplit[2].startsWith("polymorphism"))
			{
				prediction = MutationTasterClassification.polymorphism;
			}
			else
			{
				throw new Exception("bad data on line " + line);
			}
			
			muttasterFile.put(chrPosRefAlt, prediction);
		}
		s.close();
	}
	
	public Judgment classifyVariantUsingMutationTaster2Results(String chr, String pos, String ref, String alt) throws Exception
	{
		String key = chr+":"+pos+","+ref+","+alt;
		if(muttasterFile.containsKey(key))
		{
			MutationTasterClassification pc = muttasterFile.get(key);
			if(pc.equals(MutationTasterClassification.polymorphism))
			{
				return new Judgment(Classification.Benign, Method.calibrated, key, "MutationTaster2 result 'polymorphism'");
			}
			else if(pc.equals(MutationTasterClassification.disease_causing))
			{
				return new Judgment(Classification.Pathogn, Method.calibrated, key, "MutationTaster2 result 'disease_causing'");
			}
			else
			{
				throw new VariantClassificationException("Unknown MutationTaster2 result: " + pc);
			}
		}
		else
		{
			throw new VariantClassificationException("No MutationTaster2 result");
		}
	}
	
	
}
