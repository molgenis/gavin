package org.molgenis.calibratecadd.support;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.HashMap;
import java.util.Scanner;

public class LoadCADDWebserviceOutput
{

	/**
	 * "chr_pos_ref_alt" to CADD PHRED score
	 * @param caddFile
	 * @return
	 * @throws FileNotFoundException
	 */
	public static HashMap<String, Double> load(File caddFile) throws FileNotFoundException
	{
		Scanner cadd = new Scanner(caddFile);
		
		HashMap<String, Double> caddScores = new HashMap<String, Double>();
		
		String line = null;
		while(cadd.hasNextLine())
		{
			line = cadd.nextLine();
			if(line.startsWith("#"))
			{
				continue;
			}
			String[] split = line.split("\t", -1);
			caddScores.put(split[0] + "_" + split[1] + "_" + split[2] + "_" + split[3], Double.parseDouble(split[5]));
		}
		cadd.close();
		return caddScores;
	}
	
}
