package org.molgenis.calibratecadd.support;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.HashMap;
import java.util.Scanner;

import org.molgenis.data.annotation.joeri282exomes.Judgment;
import org.molgenis.data.annotation.joeri282exomes.Judgment.Classification;
import org.molgenis.data.annotation.joeri282exomes.Judgment.Method;

public class PolyPhen2Results
{
	
	enum PolyPhenClassification {
		Benign, Damaging
	}

	//chr:pos,ref,alt to PolyPhenClassification
	HashMap<String, String> polyphenFile;
	
	/**
	 * 
	 * Example:
	 * 
	 * 10	112572605	G	T	benign
	 * 14	23857459	G	A	probably damaging
	 * 1	237821276	T	C	possibly damaging
	 * 
	 * @throws FileNotFoundException 
	 * 
	 */
	public PolyPhen2Results(File polyphen2results) throws FileNotFoundException
	{
		polyphenFile = new HashMap<String, String>();
		Scanner s = new Scanner(polyphen2results);
	//	s.nextLine(); no header!
		while(s.hasNextLine())
		{
			String line = s.nextLine();
			String[] linesplit = line.split("\t", -1);
			polyphenFile.put(linesplit[0] + "_" + linesplit[1] + "_" + linesplit[2] + "_" + linesplit[3], linesplit[4]);
		}
		s.close();
	}
	
	public Judgment classifyVariantUsingPolyPhen2Results(String chr, String pos, String ref, String alt) throws Exception
	{
		String key = chr+"_"+pos+"_"+ref+"_"+alt;
		if(polyphenFile.containsKey(key))
		{
			String pc = polyphenFile.get(key);
			
			if(pc.contains("benign"))
			{
				return new Judgment(Classification.Benign, Method.calibrated, "PolyPhen2 result 'benign'");
			}
			else if(pc.contains("damaging"))
			{
				return new Judgment(Classification.Pathogn, Method.calibrated, "PolyPhen2 result 'damaging'");
			}
			else
			{
				return new Judgment(Classification.VOUS, Method.calibrated, "Unknown PolyPhen2 result: " + pc);
			}
		}
		else
		{
			return new Judgment(Classification.VOUS, Method.calibrated, "No PolyPhen2 result");
		}
	}
	
	
}
