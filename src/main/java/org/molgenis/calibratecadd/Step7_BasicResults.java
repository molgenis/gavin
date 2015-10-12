package org.molgenis.calibratecadd;

import java.io.File;
import java.io.PrintWriter;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Scanner;

import org.apache.commons.math3.stat.descriptive.moment.Mean;
import org.apache.commons.math3.stat.descriptive.rank.Median;
import org.apache.commons.math3.stat.inference.MannWhitneyUTest;

import net.mintern.primitive.Primitive;

/**
 * 
 * Read:
 * 
 * gene	chr	pos	ref	alt	group	cadd
 * IFT172	2	27680627	A	T	PATHOGENIC	28.0
 * IFT172	2	27700177	A	T	PATHOGENIC	15.99
 * IFT172	2	27693963	C	T	PATHOGENIC	36
 * 
 * Write:
 * 
 * gene	nPath	nPopul	medianPatho	medianPopul	medianDiff
 * IFT172	10	14	22.35	25.69	3.34
 *
 *(more? highest, lowest, averages...)
 *
 *
 * Example:
 * E:\Data\clinvarcadd\clinvar.patho.fix.snpeff.exac.withcadd.tsv
 * E:\Data\clinvarcadd\clinvar.patho.fix.snpeff.exac.genesumm.tsv
 *
 */
public class Step7_BasicResults {

	public static void main(String[] args) throws Exception {
		System.out.println("starting..");
		
		Scanner res = new Scanner(new File(args[0]));
		PrintWriter pw = new PrintWriter(new File(args[1]));
		
		pw.println("gene" + "\t" + "nPath" + "\t" + "nPopul" + "\t" + "medianPatho" + "\t" + "medianPopul" + "\t" + "medianDiff" + "\t" + "meanPatho" + "\t" + "meanPopul" + "\t" + "meanDiff" + "\t" + "uTestPval" + "\t" + "caddPercConfThreshold" + "\t" + "nrOfPopulAbove95PercThreshold" + "\t" + "nrOfPathoAbove95PercThreshold");
		
		NumberFormat f = new DecimalFormat("#0.00");     
		
		HashMap<String, ArrayList<String>> geneToLines = new HashMap<String, ArrayList<String>>();
		res.nextLine();
		String line = null;
		while(res.hasNextLine())
		{
			line = res.nextLine();
			String gene = line.split("\t", -1)[0];
			if(geneToLines.containsKey(gene))
			{
				geneToLines.get(gene).add(line);
			}
			else
			{
				ArrayList<String> lines = new ArrayList<String>();
				lines.add(line);
				geneToLines.put(gene, lines);
			}
		}
		
		int nrOfGenesPathGtPopPval_5perc = 0;
		int nrOfGenesPathGtPopPval_1perc = 0;
		
		for(String gene : geneToLines.keySet())
		{
			ArrayList<Double> caddPatho = new ArrayList<Double>();
			ArrayList<Double> caddPopul = new ArrayList<Double>();
			for(String lineForGene : geneToLines.get(gene))
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
					res.close();
					throw new Exception("unknown group " + group);
				}
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
			
			Median median = new Median();
			double pathoMedian = median.evaluate(caddPathoPrim);
			double populMedian = median.evaluate(caddPopulPrim);
			double medianDiff = pathoMedian-populMedian;
			
			Mean mean = new Mean();
			double pathoMean = mean.evaluate(caddPathoPrim);
			double populMean = mean.evaluate(caddPopulPrim);
			double meanDiff = pathoMean-populMean;
			
			MannWhitneyUTest utest = new MannWhitneyUTest();
			double pval = utest.mannWhitneyUTest(caddPathoPrim, caddPopulPrim);
			
			/**
			 * Get 95% confidence of pathogenicity CADD thresholds, as best we can
			 */
			double cadd95PercConfThreshold = -1;
			int nrOfPopAbove95PercThreshold = -1;
			int nrOfPathAbove95PercThreshold = -1;
			Primitive.sort(caddPathoPrim, (d1, d2) -> Double.compare(d1, d2), false);
			
			for(int i = 0; i < caddPathoPrim.length; i++)
			{
				double threshold = caddPathoPrim[i];
				int nrOfPopAboveThreshold = 0;
				for(int n = 0; n < caddPopulPrim.length; n++)
				{
					if(caddPopulPrim[n] > threshold)
					{
						nrOfPopAboveThreshold++;
					}
				}
				
				int nrOfPathAboveThreshold = caddPathoPrim.length-i;
				double ratio = (double)nrOfPopAboveThreshold/(double)nrOfPathAboveThreshold;
		
				//happens when either we truly have 95% confidence (a 1:20 ratio), or we run out of population variants
				//while there are still 1+ pathogenic variants left
				//FIXME: not quite okay yet, some artifacts due to low population variants in e.g. MECP2
				if(ratio <= 0.05)
				{
					cadd95PercConfThreshold = threshold;
					nrOfPopAbove95PercThreshold = nrOfPopAboveThreshold;
					nrOfPathAbove95PercThreshold = nrOfPathAboveThreshold;
					break;
				}
			}
			
			//if we have no threshold (population variants have full overlap with pathogenic, and not enough samples to get 1:20 ratio)
			//we assign the maximum value of the population.. not ideal but better than nothing. Just show this in the final table (counts of 0).
			if(cadd95PercConfThreshold == -1)
			{
				cadd95PercConfThreshold = Collections.max(caddPopul);
				nrOfPopAbove95PercThreshold = 0;
				nrOfPathAbove95PercThreshold = 0;
			}
			
			//to show some stats in the sysout
			if(pval < 0.05 && pathoMean > populMean)
			{
				nrOfGenesPathGtPopPval_5perc ++;
				if(pval < 0.01)
				{
					nrOfGenesPathGtPopPval_1perc++;
				}
			}
			
			//write table
			pw.println(gene + "\t" + caddPathoPrim.length + "\t" + caddPopulPrim.length + "\t" + f.format(pathoMedian) + "\t" + f.format(populMedian) + "\t" + f.format(medianDiff) + "\t" + f.format(pathoMean) + "\t" + f.format(populMean) + "\t" + f.format(meanDiff) + "\t" + pval + "\t" + cadd95PercConfThreshold + "\t" + nrOfPopAbove95PercThreshold + "\t" + nrOfPathAbove95PercThreshold);
		
		}
		
		System.out.println("total nr of genes: " + geneToLines.keySet().size());
		System.out.println("nr of genes where patho > pop, pval < 0.05: " + nrOfGenesPathGtPopPval_5perc);
		System.out.println("nr of genes where patho > pop, pval < 0.01: " + nrOfGenesPathGtPopPval_1perc);
		
		pw.flush();
		pw.close();
		
		res.close();

	}

}
