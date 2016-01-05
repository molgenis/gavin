package org.molgenis.calibratecadd;

import java.io.File;
import java.io.PrintWriter;
import java.util.HashMap;
import java.util.Scanner;

import org.molgenis.calibratecadd.support.LoadCADDWebserviceOutput;

/**
 * Example usage:
 * E:\Data\clinvarcadd\cadd_output.txt
 * E:\Data\clinvarcadd\clinvar.patho.fix.snpeff.exac.vcf.info
 * E:\Data\clinvarcadd\clinvar.patho.fix.snpeff.exac.withcadd.tsv
 * 
 * 
 * Takes the CADD output, which looks like this:
 * 
 * ## CADD v1.3 (c) University of Washington and Hudson-Alpha Institute for Biotechnology 2013-2015. All rights reserved.
 * #CHROM	POS	REF	ALT	RawScore	PHRED
 * 1	976962	C	T	6.039163	28.0
 * 1	976963	A	G	-2.588522	0.001
 * 1	976967	G	A	0.552271	7.819
 * etc
 *
 *
 * Plus all accompanying information outputted directly by Step 4:
 * 
 * gene	chr	pos	ref	alt	group
 * IFT172	2	27680627	A	T	PATHOGENIC
 * IFT172	2	27700177	A	T	PATHOGENIC
 * IFT172	2	27693963	C	T	PATHOGENIC
 * etc
 * 
 * 
 * Combine this into:
 * 
 * gene	chr	pos	ref	alt	group	cadd
 * IFT172	2	27680627	A	T	PATHOGENIC	28.0
 * IFT172	2	27700177	A	T	PATHOGENIC	15.99
 * IFT172	2	27693963	C	T	PATHOGENIC	36
 * etc
 * 
 *
 */
public class Step6_PrepForAnalysis
{

	public static void main(String[] args) throws Exception
	{
//		Scanner cadd = new Scanner(new File(args[0]));
		Scanner info = new Scanner(new File(args[1]));
		PrintWriter pw = new PrintWriter(new File(args[2]));
		
		HashMap<String, Double> caddScores = LoadCADDWebserviceOutput.load(new File(args[0]));
		
//		HashMap<String, Double> caddScores = new HashMap<String, Double>();
//		
		String line = null;
//		while(cadd.hasNextLine())
//		{
//			line = cadd.nextLine();
//			if(line.startsWith("#"))
//			{
//				continue;
//			}
//			String[] split = line.split("\t", -1);
//			caddScores.put(split[0] + "_" + split[1] + "_" + split[2] + "_" + split[3], Double.parseDouble(split[5]));
//		}
//		cadd.close();
		
		//write header of output
		pw.println("gene" + "\t" + "chr" + "\t" + "pos" + "\t" + "ref" + "\t" + "alt" + "\t" + "group" + "\t" + "cadd");
		
		//skip header of input
		info.nextLine();
		
		while(info.hasNextLine())
		{
			line = info.nextLine();
			String[] split = line.split("\t", -1);
			String key = split[1] + "_" + split[2] + "_" + split[3] + "_" + split[4];
			if(caddScores.containsKey(key))
			{
				//FIXME: need to replace '/' to prevent problems in R later on when writing plots based on gene names..
				pw.println(line.replace("/", "_") + "\t" + caddScores.get(key));
			}
			else
			{
				System.out.println("WARNING: could not get CADD score for " + key);
			}
		}
		
		info.close();
		pw.flush();
		pw.close();
		
		
	}
}
