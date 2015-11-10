package org.molgenis.calibratecadd;

import java.io.File;
import java.io.PrintWriter;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Scanner;

import org.apache.commons.lang3.StringUtils;
import org.apache.commons.math3.stat.descriptive.moment.Mean;
import org.apache.commons.math3.stat.descriptive.rank.Percentile;
import org.apache.commons.math3.stat.descriptive.rank.Percentile.EstimationType;
import org.apache.commons.math3.stat.inference.MannWhitneyUTest;

/**
 * 
 * Read:
 * (e.g. clinvar.patho.fix.snpeff.exac.genes.tsv, output of step 4)
 * 
 * gene	category	info
 * NUP107	N1	less than 2 clinvar variants available
 * IFT172	Cx	ready for cadd calibration with 96 variants
 * MT-TV	N2	0 exac variants in range MT:1606-1644
 * NDST1	T1	no exac variants below MAF 7.000599999999998E-6 known impacts: ImpactRatios [high=0.0, moderate=100.0, low=0.0, modifier=0.0]
 * 
 * 
 * And:
 * (e.g. clinvar.patho.fix.snpeff.exac.withcadd.tsv, output of step 6)
 * 
 * gene	chr	pos	ref	alt	group	cadd
 * IFT172	2	27680627	A	T	PATHOGENIC	28.0
 * IFT172	2	27700177	A	T	PATHOGENIC	15.99
 * IFT172	2	27693963	C	T	PATHOGENIC	36
 * 
 * 
 * Write out:
 * (e.g. clinvar.patho.fix.snpeff.exac.genesumm.tsv)
 * 
 * gene	nPath	nPopul	medianPatho	medianPopul	medianDiff
 * IFT172	10	14	22.35	25.69	3.34
 *
 *
 *
 * Example:
 * E:\Data\clinvarcadd\clinvar.patho.fix.snpeff.exac.genes.tsv
 * E:\Data\clinvarcadd\clinvar.patho.fix.snpeff.exac.withcadd.tsv
 * E:\Data\clinvarcadd\clinvar.patho.fix.snpeff.exac.genesumm.tsv
 *
 */
public class Step7_BasicResults {

	public static void main(String[] args) throws Exception {
		System.out.println("starting..");
		
		/**
		 * read gene info and put in map
		 */
		HashMap<String, String> geneToInfo = new HashMap<String, String>();
		Scanner geneInfoScanner = new Scanner(new File(args[0]));
		geneInfoScanner.nextLine(); //skip header
		String line = null;
		while(geneInfoScanner.hasNextLine())
		{
			line = geneInfoScanner.nextLine();
			String[] split = line.split("\t", -1);
			line = line.substring(line.indexOf('\t')+1, line.length());
			geneToInfo.put(split[0], line);
		}
		geneInfoScanner.close();
		
		/**
		 * read variant + cadd data and put in map
		 */
		Scanner variantsWithCADDScanner = new Scanner(new File(args[1]));
		HashMap<String, ArrayList<String>> geneToVariantAndCADD = new HashMap<String, ArrayList<String>>();
		variantsWithCADDScanner.nextLine();
		line = null;
		while(variantsWithCADDScanner.hasNextLine())
		{
			line = variantsWithCADDScanner.nextLine();
			String gene = line.split("\t", -1)[0];
			if(geneToVariantAndCADD.containsKey(gene))
			{
				geneToVariantAndCADD.get(gene).add(line);
			}
			else
			{
				ArrayList<String> lines = new ArrayList<String>();
				lines.add(line);
				geneToVariantAndCADD.put(gene, lines);
			}
		}
		variantsWithCADDScanner.close();
		
		/**
		 * process everything and write out
		 */
		PrintWriter pw = new PrintWriter(new File(args[2]));
		pw.println("Gene" + "\t" + "Category" + "\t" + "Chr" + "\t" + "Start" + "\t" + "End" + "\t" + "NrOfPopulationVariants" + "\t" + "NrOfPathogenicVariants" + "\t" + "PathoMAFThreshold" + "\t" + "PopImpactHighPerc" + "\t" + "PopImpactModeratePerc" + "\t" + "PopImpactLowPerc" + "\t" + "PopImpactModifierPerc" + "\t" + "PathoImpactHighPerc" + "\t" + "PathoImpactModeratePerc" + "\t" + "PathoImpactLowPerc" + "\t" + "PathoImpactModifierPerc" + "\t" + "PopImpactHighEq" + "\t" + "PopImpactModerateEq" + "\t" + "PopImpactLowEq" + "\t" + "PopImpactModifierEq" + "\t" + "NrOfCADDScoredPopulationVars" + "\t" + "NrOfCADDScoredPathogenicVars" + "\t" + "MeanPopulationCADDScore" + "\t" + "MeanPathogenicCADDScore" + "\t" + "MeanDifference" + "\t" + "UTestPvalue" + "\t" + "Sens95thPerCADDThreshold" +"\t" + "Spec95thPerCADDThreshold");
		NumberFormat f = new DecimalFormat("#0.00");     
		
		int nrOfGenesPathGtPopPval_5perc = 0;
		int nrOfGenesPathGtPopPval_1perc = 0;
		
		for(String gene : geneToInfo.keySet())
		{
			if(!geneToVariantAndCADD.containsKey(gene))
			{
				pw.println(gene + "\t" + geneToInfo.get(gene) + StringUtils.repeat("\t" + Step4_MatchingVariantsFromExAC.NA, 8));
			}
			else
			{
				ArrayList<Double> caddPatho = new ArrayList<Double>();
				ArrayList<Double> caddPopul = new ArrayList<Double>();
				for(String lineForGene : geneToVariantAndCADD.get(gene))
				{
					String[] split = lineForGene.split("\t", -1);
					String group = split[5];
					double cadd = Double.parseDouble(split[6]);
					if(group.equals("PATHOGENIC"))
					{
						caddPatho.add(cadd);
					}
					else if(group.equals("POPULATION"))
					{
						caddPopul.add(cadd);
					}
					else
					{
						pw.close();
						throw new Exception("unknown group " + group);
					}
				}
				
				//it can happen that variants for one group did not pass CADD webservice, e.g. for PRRT2 we have only 1 population variant and when that fails, we have cannot assess...
				//replace 'Cx' with 'N3'
				if(caddPatho.size() == 0 || caddPopul.size() == 0)
				{
					pw.println(gene + "\t" + "N3" + geneToInfo.get(gene).substring(2, geneToInfo.get(gene).length()) + "\t" + caddPopul.size() + "\t" + caddPatho.size() + StringUtils.repeat("\t" + Step4_MatchingVariantsFromExAC.NA, 6));
					continue;
				}
				
				double[] caddPathoPrim = new double[caddPatho.size()];
				for(int i = 0; i < caddPatho.size(); i++)
				{
					caddPathoPrim[i] = caddPatho.get(i);
				}
				
				double[] caddPopulPrim = new double[caddPopul.size()];
				for(int i = 0; i < caddPopul.size(); i++)
				{
					caddPopulPrim[i] = caddPopul.get(i);
				}
				
				Mean mean = new Mean();
				double pathoMean = mean.evaluate(caddPathoPrim);
				double populMean = mean.evaluate(caddPopulPrim);
				double meanDiff = pathoMean-populMean;
				
				MannWhitneyUTest utest = new MannWhitneyUTest();
				double pval = utest.mannWhitneyUTest(caddPathoPrim, caddPopulPrim);
				
				//get thresholds for 95% sensitivity and 95% specificity
				//sensitive: we catch 95% of the known pathogenic variants, no matter how many population variants we also include when using this threshold
				//specific: we only allow a 'top 5%' error (=finding population variants) in finding pathogenic variants, no matter how many pathogenic variants we have above this threshold right now
				Percentile perc = new Percentile().withEstimationType(EstimationType.R_7);
				double sensThres = perc.evaluate(caddPathoPrim, 5);
				double specThres = perc.evaluate(caddPopulPrim, 95);
				
				String cat = null;
				//to show some stats in the sysout
				if(pval <= 0.05 && pathoMean > populMean)
				{
					cat = "C2";
					nrOfGenesPathGtPopPval_5perc ++;
					if(pval <= 0.01)
					{
						cat = "C1";
						nrOfGenesPathGtPopPval_1perc++;
					}
				}
				
				if(pval > 0.05)
				{
					if(caddPathoPrim.length < 5 || caddPopulPrim.length < 5)
					{
						cat = "C3";
					}
					else
					{
						cat = "C4";
					}
				}
				
				
				if(cat == null)
				{
					cat = "C5";
				}
				
				//write table
				pw.println(gene + "\t" + cat + geneToInfo.get(gene).substring(2, geneToInfo.get(gene).length()) + "\t" + caddPopulPrim.length + "\t" + caddPathoPrim.length + "\t" + f.format(populMean) + "\t" + f.format(pathoMean) + "\t" + f.format(meanDiff) + "\t" + pval + "\t" + f.format(sensThres) + "\t" + f.format(specThres));
			}
		}
		
		System.out.println("total nr of genes: " + geneToVariantAndCADD.keySet().size());
		System.out.println("nr of genes where patho > pop, pval < 0.05: " + nrOfGenesPathGtPopPval_5perc);
		System.out.println("nr of genes where patho > pop, pval < 0.01: " + nrOfGenesPathGtPopPval_1perc);
		
		pw.flush();
		pw.close();

	}

}
