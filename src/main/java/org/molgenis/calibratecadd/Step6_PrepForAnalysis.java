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
 * IFT172	2	27680627	A	T	PATHOGENIC	splice_region_variant&intron_variant	LOW
 * IFT172	2	27700177	A	T	PATHOGENIC	missense_variant	MODERATE
 * IFT172	2	27693963	C	T	PATHOGENIC	splice_acceptor_variant&intron_variant	HIGH
 * IFT172	2	27667971	TCCTGTG	T	POPULATION	frameshift_variant&splice_acceptor_variant&splice_region_variant&intron_variant	HIGH	36.0
 * etc
 * 
 * 
 * Combine / clean this into:
 * 
 * gene	chr	pos	ref	alt	group	cadd
 * IFT172	2	27680627	A	T	PATHOGENIC	28.0	splice_region_variant	LOW
 * IFT172	2	27700177	A	T	PATHOGENIC	15.99	missense_variant	MODERATE
 * IFT172	2	27693963	C	T	PATHOGENIC	36	splice_acceptor_variant	HIGH
 * IFT172	2	27667971	TCCTGTG	T	POPULATION	frameshift_variant	HIGH	36.0
 * etc
 * 
 *
 */
public class Step6_PrepForAnalysis
{

	public static void main(String[] args) throws Exception
	{
		HashMap<String, Double> caddScores = LoadCADDWebserviceOutput.load(new File(args[0]));
		Scanner variants = new Scanner(new File(args[1]));
		PrintWriter pw = new PrintWriter(new File(args[2]));
		
		String line = null;
		
		//write header of output
		pw.println("gene" + "\t" + "chr" + "\t" + "pos" + "\t" + "ref" + "\t" + "alt" + "\t" + "group" + "\t" + "effect" + "\t" + "impact" + "\t" + "cadd");
		
		//skip header of input
		variants.nextLine();
		
		while(variants.hasNextLine())
		{
			line = variants.nextLine();
			String[] split = line.split("\t", -1);
			String gene = split[0];
			String chr = split[1];
			String pos = split[2];
			String ref = split[3];
			String alt = split[4];
			String group = split[5];
			String effect = split[6];
			String impact = split[7];

			effect = effect.substring(0, (effect.indexOf("&") == -1 ? effect.length() : effect.indexOf("&")));
			String printMe = gene + "\t" + chr + "\t" + pos + "\t" + ref + "\t" + alt + "\t" + group + "\t" + effect + "\t" + impact;

			String key = chr + "_" + pos + "_" + ref + "_" + alt;
			if(caddScores.containsKey(key))
			{
				//FIXME: need to replace '/' to prevent problems in R later on when writing plots based on gene names..
				pw.println(printMe.replace("/", "_") + "\t" + caddScores.get(key));
			}
			else
			{
				String trimmedRefAlt = LoadCADDWebserviceOutput.trimRefAlt(ref, alt, "_");
				key = chr + "_" + pos + "_" + trimmedRefAlt;
				if(caddScores.containsKey(key))
				{
					System.out.println("RESOLVED by trimming ref alt " + ref + "_" + alt + " to " + trimmedRefAlt);
					pw.println(printMe.replace("/", "_") + "\t" + caddScores.get(key));
				}
				else
				{
					System.out.println("WARNING: could not get CADD score for " + key);
				}


			}
		}

		variants.close();
		pw.flush();
		pw.close();
		
		
	}
}
