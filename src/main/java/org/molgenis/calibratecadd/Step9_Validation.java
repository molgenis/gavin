package org.molgenis.calibratecadd;

import java.io.File;
import java.io.FileNotFoundException;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.HashMap;
import java.util.Scanner;

public class Step9_Validation
{

	/**
	 * Takes a Cartagenia exported file with e.g. "likely benign" variants as args[0].
	 * Example lines:
	 * 
	 * chromosome	start	stop	ref	alt	variant_type	location	effect	gene	transcript	exon	c_nomen	p_nomen	dbsnp
	 * 1	78383891	78383891	G	A	snp	exonic	nonsynonymous	NEXN	NM_144573.3	5	c.380G>A	p.R127H	NULL
	 * 1	78395131	78395131	A	C	snp	exonic	nonsynonymous	NEXN	NM_144573.3	9	c.995A>C	p.E332A	rs201763096
	 * 
	 * 
	 * 
	 * Takes a file with CADD score results from the webservice as args[1].
	 * Example lines:
	 * 
	 * ## CADD v1.3 (c) University of Washington and Hudson-Alpha Institute for Biotechnology 2013-2015. All rights reserved.
	 * #CHROM	POS	REF	ALT	RawScore	PHRED
	 * 1	78383891	G	A	6.630134	32
	 * 1	78392121	T	C	1.153591	11.50
	 * 
	 * 
	 * 
	 * Takes the 'gene summary' output of Step 7 as args[2].
	 * Example lines:
	 * 
	 * Gene    NrOfPathoVars   NrOfPopulVars   MeanPathoCADD   MeanPopulCADD   MeanDiff        UTestPvalue     Sens95perCADDthresh     Spec95perCADDthresh
	 * IFT172  15      82      28.15   26.57   1.58    0.3017541066910733      0.00    41.00
	 * SPINT2  6       7       24.85   23.53   1.32    0.6682351417952495      22.20   32.00
	 * NDST1   4       11      29.57   27.15   2.42    0.5138908912577136      24.20   34.00
	 * 
	 * 
	 * 
	 * Output, example:
	 * 85% of variants was found above the 95% sensitivity threshold
	 * 12% of variants was found above the 95% specificity threshold
	 *
	 * @param args
	 * @throws FileNotFoundException
	 */
	public static void main(String[] args) throws FileNotFoundException
	{
		//1_78383891_G_A -> 32
		HashMap<String, Double> caddScores = new HashMap<String, Double>();
		Scanner s = new Scanner(new File(args[1]));
		String line = null;
		while(s.hasNextLine())
		{
			line = s.nextLine();
			if(line.startsWith("#"))
			{
				continue;
			}
			String[] split = line.split("\t", -1);
			caddScores.put(split[0]+"_"+split[1]+"_"+split[2]+"_"+split[3], Double.parseDouble(split[5]));
		}
		s.close();
		
		//NDST1 -> 24.20
		//NDST1 -> 34.00
		HashMap<String, Double> geneThrSens = new HashMap<String, Double>();
		HashMap<String, Double> geneThrSpec = new HashMap<String, Double>();
		s = new Scanner(new File(args[2]));
		line = null;
		while(s.hasNextLine())
		{
			line = s.nextLine();
			if(line.startsWith("Gene"))
			{
				continue;
			}
			String[] split = line.split("\t", -1);
			geneThrSens.put(split[0], Double.parseDouble(split[7]));
			geneThrSpec.put(split[0], Double.parseDouble(split[8]));
		}
		s.close();
		
		//now check cartagenia MVL export
		s = new Scanner(new File(args[0]));
		line = null;
		int nrOfVariantsOverSensThr = 0;
		int nrOfVariantsOverSpecThr = 0;
		int nrOfVariantsEvaluated = 0;
		int nrOfVariantsTotal = 0;
		while(s.hasNextLine())
		{
			line = s.nextLine();
			String[] split = line.split("\t", -1);
			String key = split[0]+"_"+split[1]+"_"+split[3]+"_"+split[4];
			String gene = split[8];
			nrOfVariantsTotal++;
			
			if(caddScores.containsKey(key) && geneThrSens.containsKey(gene))
			{
				//just a check...
				if(geneThrSens.get(gene) > geneThrSpec.get(gene))
				{
					System.out.println("WARNING: sensitiviy threshold higher than specificity threshold for " + gene + " !!!");
				}
				
				nrOfVariantsEvaluated++;
				
				double score = caddScores.get(key);
				if(score > geneThrSens.get(gene))
				{
					nrOfVariantsOverSensThr++;
				}
				if(score > geneThrSpec.get(gene))
				{
					nrOfVariantsOverSpecThr++;
				}
			}
			
		}
		
		System.out.println("cartagenia file: " + args[0]);
		System.out.println("number of variants total: " + nrOfVariantsTotal);
		System.out.println("number of variants evaluated: " + nrOfVariantsEvaluated);
		System.out.println("number of variants over sensitivity threshold: " + nrOfVariantsOverSensThr);
		System.out.println("number of variants over specificity threshold: " + nrOfVariantsOverSpecThr);
		
		NumberFormat f = new DecimalFormat("#0.00");   
		System.out.println("percentage of evaluated variants over sensitivity threshold: " + f.format(((double)nrOfVariantsOverSensThr/(double)nrOfVariantsEvaluated)*100.0));
		System.out.println("percentage of evaluated variants over specificity threshold: " + f.format(((double)nrOfVariantsOverSpecThr/(double)nrOfVariantsEvaluated)*100.0));
		
		s.close();

	}

}
