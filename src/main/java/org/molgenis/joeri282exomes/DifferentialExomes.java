package org.molgenis.joeri282exomes;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Scanner;

import org.apache.commons.math3.stat.descriptive.moment.Mean;
import org.apache.commons.math3.stat.inference.TestUtils;
import org.molgenis.data.Entity;
import org.molgenis.data.annotator.tabix.TabixVcfRepository;
import org.molgenis.data.vcf.VcfRepository;

public class DifferentialExomes
{
	private File vcfFile;
	private File patientGroups;
	private TabixVcfRepository vcfReader;
	private HashMap<String, String> sampleToGroup;

	public DifferentialExomes(File vcfFile, File patientGroups) throws IOException
	{

		// TODO: checks and stuff

		this.vcfFile = vcfFile;
		this.patientGroups = patientGroups;
		this.vcfReader = new TabixVcfRepository(vcfFile, "patients");

		HashMap<String, String> sampleToGroup = new HashMap<String, String>();
		Scanner s = new Scanner(patientGroups);
		String line = null;
		while (s.hasNextLine())
		{
			line = s.nextLine();
			String[] lineSplit = line.split("\t", -1);
			sampleToGroup.put(lineSplit[0], lineSplit[1]);
		}
		s.close();

		this.sampleToGroup = sampleToGroup;
	}

	public void start() throws Exception
	{

		PrintWriter pw = new PrintWriter("/Users/jvelde/Desktop/joeri282exomes/DiffExOut.tsv", "UTF-8");
		HashMap<String, HashMap<String, Boolean>> geneToSampleToLOF = getSampleLOF();

		for (String sample : sampleToGroup.keySet())
		{
			pw.print("\t" + sample);
		}
		pw.print("\n");
		pw.flush();

		for (String gene : geneToSampleToLOF.keySet())
		{
			StringBuilder sb = new StringBuilder();
			sb.append(gene);

			for (String sample : sampleToGroup.keySet())
			{
				sb.append("\t" + geneToSampleToLOF.get(gene).get(sample));
			}

			System.out.println("GENE: " + gene);
			// calculateGroupDiff(geneToSampleToLOF.get(gene));
			calculateGroupChiSq(geneToSampleToLOF.get(gene));
			System.out.println("----");

			pw.println(sb.toString());
			pw.flush();
		}

		pw.close();

	}

	private void calculateGroupChiSq(HashMap<String, Boolean> sampleToLOF)
	{
		HashMap<String, ArrayList<Boolean>> valuesPerGroup = new HashMap<String, ArrayList<Boolean>>();

		// example from: http://commons.apache.org/proper/commons-math/userguide/stat.html
		// long[] observed = {10, 9, 11};
		// double[] expected = {10.1, 9.8, 10.3};
		// System.out.println(TestUtils.chiSquareTest(expected, observed));
		//
		// in our case:
		// observed = { nr of LoF in patients, nr of no-LoF in patients};
		// expected = { nr of LoF in controls, nr of no-LoF in controls};

		// in R:
		// LOF <- matrix(c(10, 20, 100, 150),nrow = 2,dimnames =list(c("LoF", "No LoF"),c("Affected", "Control")))
		// fisher.test(LOF, alternative="greater")
		// p-value = 0.8164 :-)
		//
		// java:
		// double pval = FishersExactTest.run(10, 20, 100, 150);
		// pval 0.8163955241091263
		System.out.println(pval);

	}

	/**
	 * Check differences between groups for 1 gene
	 * 
	 * @param sampleToLOF
	 */
	private void calculateGroupDiff(HashMap<String, Boolean> sampleToLOF)
	{
		HashMap<String, ArrayList<Boolean>> valuesPerGroup = new HashMap<String, ArrayList<Boolean>>();

		double[] allValues = new double[sampleToLOF.size()];
		int i = 0;
		for (String sample : sampleToGroup.keySet())
		{
			allValues[i++] = sampleToLOF.get(sample) == true ? 1.0 : 0.0;
			String group = sampleToGroup.get(sample);

			// put samples into different buckets depending on their group 'annotation'
			if (valuesPerGroup.containsKey(group))
			{
				valuesPerGroup.get(group).add(sampleToLOF.get(sample));
			}
			else
			{
				ArrayList<Boolean> values = new ArrayList<Boolean>();
				values.add(sampleToLOF.get(sample));
				valuesPerGroup.put(group, values);
			}
		}

		// since we have seen all values, why not keep them to calculate the mean:
		Mean m = new Mean();
		double mean = m.evaluate(allValues);

		System.out.println("mean for all samples is " + mean);

		// now calculate mean for each of the buckets
		i = 0;
		for (String group : valuesPerGroup.keySet())
		{
			double[] groupValues = new double[valuesPerGroup.get(group).size()];
			for (Boolean value : valuesPerGroup.get(group))
			{
				groupValues[i++] = value == true ? 1.0 : 0.0;
			}
			i = 0;
			m = new Mean();
			double groupMean = m.evaluate(groupValues);

			System.out.println("mean for group " + group + " is " + groupMean);

		}

	}

	public HashMap<String, HashMap<String, Boolean>> getSampleLOF() throws Exception
	{
		HashMap<String, HashMap<String, Boolean>> geneToSampleToLOF = new HashMap<String, HashMap<String, Boolean>>();

		Iterator<Entity> vcfIter = vcfReader.iterator();

		// keep track of the 'line' in which variant was seen
		int index = 0;

		vcfIterNext: while (vcfIter.hasNext())
		{
			index++;

			if (index % 1000 == 0)
			{
				System.out.println("Seen " + index + " variants..");
			}

			if (index % 10000 == 0)
			{
				System.out.println("DEV !! quitting...");
				break;
			}

			Entity record = vcfIter.next();

			// System.out.println("looking at " + record.toString());

			String chr = record.getString("#CHROM");
			String pos = record.getString("POS");
			String ref = record.getString("REF");
			String alt = record.getString("ALT");

			String[] altSplit = alt.split(",", -1);

			if (altSplit.length > 1)
			{
				// TODO
				// System.out.println("WARNING: skipping multi allelic PATIENT VARIANT for now..");
				continue;
			}

			String[] multiAnn = record.getString("INFO_ANN").split(",");

			for (int i = 0; i < multiAnn.length; i++)
			{
				String ann = multiAnn[i];
				String[] annSplit = ann.split("\\|", -1);
				String impact = annSplit[2];

				if (!impact.equals("HIGH"))
				{
					continue vcfIterNext;
				}

				String gene = annSplit[3];
				if (gene.isEmpty())
				{
					throw new Exception("reminder: gene can be empty??");
				}

				Iterable<Entity> sampleEntities = record.getEntities(VcfRepository.SAMPLES);
				for (Entity sample : sampleEntities)
				{

					String genotype = sample.get("GT").toString();

					if (genotype.equals("./."))
					{
						continue;
					}

					if (genotype.equals("1/1"))
					{
						// String sampleName = sample.get("NAME").toString().split("_")[2];
						String sampleName = sample.get("ORIGINAL_NAME").toString();
						addGeneSampleLOF(geneToSampleToLOF, gene, sampleName);

					}

				}

			}

		}

		for (String gene : geneToSampleToLOF.keySet())
		{
			// System.out.println("gene " + gene);
			for (String sample : sampleToGroup.keySet())
			{
				if (!geneToSampleToLOF.get(gene).containsKey(sample))
				{
					// System.out.println("adding sample because missing " + sample);
					geneToSampleToLOF.get(gene).put(sample, false);
				}
			}
		}

		return geneToSampleToLOF;

	}

	private void addGeneSampleLOF(HashMap<String, HashMap<String, Boolean>> geneToSampleToLOF, String gene,
			String sample) throws Exception
	{
		if (geneToSampleToLOF.keySet().contains(gene))
		{
			if (geneToSampleToLOF.get(gene).containsKey(sample))
			{
				// okay - could be multiple hits per gene/sample ofc.
				// throw new Exception("Already contains " + sample);
			}
			geneToSampleToLOF.get(gene).put(sample, true);
		}
		else
		{
			HashMap<String, Boolean> sampleToLOF = new HashMap<String, Boolean>();
			sampleToLOF.put(sample, true);
			geneToSampleToLOF.put(gene, sampleToLOF);
		}
	}

}
