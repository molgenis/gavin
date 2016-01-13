package org.molgenis.calibratecadd.support;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.HashMap;
import java.util.Scanner;

import org.apache.commons.lang3.StringUtils;

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
	
	
	/**
	 * AT ATT -> A AT
	 * ATGTG ATG -> ATG A
	 * ATGTG ATGTGTGTG -> A ATGTG
	 * GATAT GAT -> GAT G
	 */
	public static String trimRefAlt(String ref, String alt, String sep)
	{
		
		char[] refRev = StringUtils.reverse(ref).toCharArray();
		char[] altRev = StringUtils.reverse(alt).toCharArray();
		
		int nrToDelete = 0;
		for(int i = 0; i < refRev.length; i++)
		{
			char refBase = refRev[i];
			char altBase = altRev[i];
			if(refBase == altBase)
			{
				nrToDelete++;
			}
			else
			{
				break;
			}
		}
		String newRef = ref.substring(0, ref.length()-nrToDelete);
		String newAlt = alt.substring(0, alt.length()-nrToDelete);
		
//		System.out.println("turning '" + ref + " " + alt + "' into '" + newRef + " " + newAlt + "'");
		return newRef + sep + newAlt;
	}
	
	public static void main(String[] args)
	{
		trimRefAlt("AT", "ATT", "\t");
		trimRefAlt("ATGTG", "ATG", "\t");
		trimRefAlt("ATGTG", "ATGTGTGTG", "\t");
		trimRefAlt("GATAT", "GAT", "\t");
		trimRefAlt("A", "T", "\t");
		trimRefAlt("AT", "TA", "\t");
		/**
		turning 'AT ATT' into 'A AT'
		turning 'ATGTG ATG' into 'ATG A'
		turning 'ATGTG ATGTGTGTG' into 'A ATGTG'
		turning 'GATAT GAT' into 'GAT G'
		turning 'A T' into 'A T'
		turning 'AT TA' into 'AT TA'
		 */
	}
	
}
