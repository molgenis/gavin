package org.molgenis.calibratecadd;

import java.io.File;
import java.io.PrintWriter;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Scanner;

import org.apache.commons.math3.stat.descriptive.moment.Mean;
import org.apache.commons.math3.stat.descriptive.rank.Median;

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
		Scanner res = new Scanner(new File(args[0]));
		PrintWriter pw = new PrintWriter(new File(args[1]));
		
		pw.println("gene" + "\t" + "nPath" + "\t" + "nPopul" + "\t" + "medianPatho" + "\t" + "medianPopul" + "\t" + "medianDiff" + "\t" + "meanPatho" + "\t" + "meanPopul" + "\t" + "meanDiff");
		
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
			
			pw.println(gene + "\t" + caddPathoPrim.length + "\t" + caddPopulPrim.length + "\t" + f.format(pathoMedian) + "\t" + f.format(populMedian) + "\t" + f.format(medianDiff) + "\t" + f.format(pathoMean) + "\t" + f.format(populMean) + "\t" + f.format(meanDiff));
		
		}
		
		pw.flush();
		pw.close();

	}

}
