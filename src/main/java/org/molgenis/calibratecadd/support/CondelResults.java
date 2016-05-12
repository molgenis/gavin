package org.molgenis.calibratecadd.support;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.HashMap;
import java.util.Scanner;

import org.molgenis.data.annotation.joeri282exomes.Judgment;
import org.molgenis.data.annotation.joeri282exomes.Judgment.Classification;
import org.molgenis.data.annotation.joeri282exomes.Judgment.Method;

public class CondelResults
{
	
	enum CondelClassification {
		D, N
	}

	//chr:pos,ref,alt to CondelClassification
	HashMap<String, CondelClassification> condelFile;
	
	/**
	 * 
	 * Example:
	 * 
	 * CHR	START	SYMBOL	REF	ALT	MA	FATHMM	CONDEL	CONDEL_LABEL
	 * 2	220285283	DES	C	G	2.075	-3.73	0.574210916504	D
	 * 18	28660261	DSC2	C	T	2.67	1.14	0.528679653055	D
	 * X	31676248	DMD	T	C	1.24	1.11	0.444478870231	N
	 * 
	 * May be conflicting, e.g.
	 * 2	179423145	TTN	C	T	2.31	0.24	0.543126551902	D
	 * 2	179423145	TTN	C	T	2.31	0.24	0.481514229992	N
	 * 
	 */
	public CondelResults(File condelResults) throws Exception
	{
		condelFile = new HashMap<String, CondelClassification>();
		Scanner s = new Scanner(condelResults);
		s.nextLine(); //skip header
		while(s.hasNextLine())
		{
			String line = s.nextLine();
			String[] linesplit = line.split("\t", -1);
			
			String chrPosRefAlt = linesplit[0] + ":" + linesplit[1] + "," + linesplit[3] + "," + linesplit[4];
			CondelClassification prediction = null;
			if(linesplit[8].equals("D"))
			{
				prediction = CondelClassification.D;
			}
			else if(linesplit[8].equals("N"))
			{
				prediction = CondelClassification.N;
			}
			else
			{
				throw new Exception("bad data on line " + line);
			}
			if(condelFile.containsKey(chrPosRefAlt) && !condelFile.get(chrPosRefAlt).equals(prediction))
			{
				System.out.println("WARNING: conflicting classification " + line + ", not adding and removing previous entry too");
				condelFile.remove(chrPosRefAlt);
			}
			else
			{
				condelFile.put(chrPosRefAlt, prediction);
			}
			
		}
		s.close();
	}
	
	public Judgment classifyVariantUsingCondelResults(String chr, String pos, String ref, String alt) throws Exception
	{
		String key = chr+":"+pos+","+ref+","+alt;
		if(condelFile.containsKey(key))
		{
			CondelClassification pc = condelFile.get(key);
			if(pc.equals(CondelClassification.N))
			{
				return new Judgment(Classification.Benign, Method.calibrated, key, "Condel result 'N'");
			}
			else if(pc.equals(CondelClassification.D))
			{
				return new Judgment(Classification.Pathogn, Method.calibrated, key, "Condel result 'D'");
			}
			else
			{
				return new Judgment(Classification.VOUS, Method.calibrated, key, "Unknown Condel result: " + pc);
			}
		}
		else
		{
			return new Judgment(Classification.VOUS, Method.calibrated, key, "No Condel result");
		}
	}
	
	
}
