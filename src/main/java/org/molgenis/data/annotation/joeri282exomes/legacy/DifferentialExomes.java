package org.molgenis.data.annotation.joeri282exomes;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.UnsupportedEncodingException;
import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Date;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.Scanner;
import java.util.Set;

import org.apache.commons.math3.stat.descriptive.moment.Mean;
import org.molgenis.caddtlmapping.mannwhitney.Helper;
import org.molgenis.data.Entity;
import org.molgenis.data.annotator.tabix.TabixVcfRepository;
import org.molgenis.data.vcf.VcfRepository;
import org.molgenis.r.ROutputHandler;
import org.molgenis.r.RScriptExecutor;
import org.molgenis.r.StringROutputHandler;

public class DifferentialExomes
{
	private File vcfFile;
	private File patientGroups;
	private File outputDataFrame;
	private TabixVcfRepository vcfReader;
	private HashMap<String, String> sampleToGroup;
	private Set<String> groups;
	
	//TODO remove
	boolean devMode = false;

	public DifferentialExomes(File vcfFile, File patientGroups) throws IOException
	{

		// TODO: checks and stuff

		this.vcfFile = vcfFile;
		this.patientGroups = patientGroups;
		this.vcfReader = new TabixVcfRepository(vcfFile, "patients");
		Set<String> groups = new HashSet<String>();

		HashMap<String, String> sampleToGroup = new HashMap<String, String>();
		Scanner s = new Scanner(patientGroups);
		String line = null;
		while (s.hasNextLine())
		{
			line = s.nextLine();
			String[] lineSplit = line.split("\t", -1);
			sampleToGroup.put(lineSplit[0], lineSplit[1]);
			groups.add(lineSplit[1]);
		}
		s.close();

		this.sampleToGroup = sampleToGroup;
		this.groups = groups;
		
		DateFormat dateFormat = new SimpleDateFormat("yyyy-MM-dd-HH-mm-ss");
		Date date = new Date();
		String inheritance = "recessive";
		String impact = "high";
		this.outputDataFrame = new File("fisherdiffexomes_" + dateFormat.format(date) + "_" + inheritance + "_"
				+ impact + ".tsv");
	}

	public void start() throws Exception
	{

		PrintWriter pw = new PrintWriter("/Users/jvelde/Desktop/joeri282exomes/DiffExOut.tsv", "UTF-8");
		
		
		LinkedHashMap<String, String> sequenceFeatureLocations = new LinkedHashMap<String, String>();
		
		
		HashMap<String, HashMap<String, Boolean>> geneToSampleToLOF = getSampleLOF(sequenceFeatureLocations);

		for (String sample : sampleToGroup.keySet())
		{
			pw.print("\t" + sample);
		}
		pw.print("\n");
		pw.flush();
		
		//gene, to group, to pval
		HashMap<String, HashMap<String, Double>> pvals = new HashMap<String, HashMap<String, Double>>();
		
		

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
			
			//do statistics per gene
			HashMap<String, Double> groupToPval = calculateGroupChiSq(geneToSampleToLOF.get(gene), gene);
			
			pvals.put(gene, groupToPval);
			
	//		System.out.println("----");

			pw.println(sb.toString());
			pw.flush();
		}

		pw.close();
		
		
		
		System.out.println("Writing output dataframe to " + outputDataFrame.getAbsolutePath());
		PrintWriter pwDF = new PrintWriter(outputDataFrame, "UTF-8");
		
		StringBuilder header = new StringBuilder();
		header.append("SNP\tCHR\tBP");
		for(String group : groups)
		{
			header.append("\t" + "P_" + group);
		}
		pwDF.println(header.toString());
		for (String sequenceFeature : pvals.keySet())
		{
			String loc = Helper.changeChromLetterToNumber((sequenceFeatureLocations.get(sequenceFeature)));
			
			StringBuilder line = new StringBuilder();
			line.append(sequenceFeature + "\t" + loc);
			for(String group : groups)
			{
				line.append("\t" + pvals.get(sequenceFeature).get(group));
			}
			pwDF.println(line.toString());
			
		}
		pwDF.flush();
		pwDF.close();

		// SNP CHR BP P
		// rs1 1 100000 0.9148
		// rs2 1 200000 0.9371
		// rs3 1 400000 0.2861
		// rs4 1 300000 0.8304
		// rs5 1 500000 0.6417
		// rs6 1 600000 0.5191

		String scriptLocation = outputDataFrame.getAbsolutePath() + ".R";
		String plotLocation = outputDataFrame.getAbsolutePath() + ".png";
		System.out.println("Creating plot at " + plotLocation + " using " + scriptLocation);
		Rscript(scriptLocation, plotLocation, outputDataFrame.getCanonicalPath());
		

	}

	private HashMap<String, Double> calculateGroupChiSq(HashMap<String, Boolean> sampleToLOF, String geneName)
	{
		HashMap<String, Double> groupToPval = new HashMap<String, Double>();
		
		
//		System.out.println("sampleToLOF size " + sampleToLOF.size());
		HashMap<String, ArrayList<Boolean>> groupValues = new HashMap<String, ArrayList<Boolean>>();
		
		for (String sample : sampleToGroup.keySet())
		{
			String group = sampleToGroup.get(sample);

			// put samples into different buckets depending on their group 'annotation'
			if (groupValues.containsKey(group))
			{
				groupValues.get(group).add(sampleToLOF.get(sample));
//				System.out.println("Adding " + sampleToLOF.get(sample) + " for sample " + sample + " in group " + group);
			}
			else
			{
				ArrayList<Boolean> values = new ArrayList<Boolean>();
				values.add(sampleToLOF.get(sample));
				groupValues.put(group, values);
//				System.out.println("Creating " + sampleToLOF.get(sample) + " for sample " + sample + " in group " + group);
			}
		}
		
		
		//group values total count?
		int totalCount = 0;
		for(String key : groupValues.keySet())
		{
			totalCount += groupValues.get(key).size();	
		}
//		System.out.println("groupValues total size " + totalCount);
		
		
		
		//test 1 group vs all the others by counting the number of TRUE and FALSE
		for(String group : groupValues.keySet())
		{
//			System.out.println("Going to test " + group);
			int a = 0;
			int b = 0;
			int c = 0;
			int d = 0;
			
			for(Boolean value : groupValues.get(group))
			{
				if(value == true)
				{
					a++;
				}
				else
				{
					b++;
				}
			}
			
//			System.out.println("counts: " + a + "x TRUE, " + b + "x FALSE");
			
			
			// now, iterate over ALL groups again, but EXCLUDE the group we're currently looking at
			
			for(String otherGroup : groupValues.keySet())
			{
//				System.out.println("Against other group: " + otherGroup);
				if(otherGroup.equals(group))
				{
//					System.out.println("GROUP ITSELF, CONTINUE");
					continue;
				}
				for(Boolean value : groupValues.get(otherGroup))
				{
					if(value == true)
					{
						c++;
					}
					else
					{
						d++;
					}
				}
//				System.out.println("othergroup counts: " + c + "x TRUE, " + d + "x FALSE");
			}
			
			// fisher exact test
			double pval = FishersExactTest.run(a, b, c, d);
			
			groupToPval.put(group, pval);
			
			double lod = -Math.log10(pval);
			
//			System.out.println("Fisher's Exact Test on " + group + " (using: " + a +", " + b +", " + c +", " + d +"),  results in pval " + pval + " (LOD: "+lod+")");

//			if(lod > 7)
//			{
//				System.out.println("lod >7 "+geneName+", interesting: " + group + " stands out with LOD " + lod);
//			}
//			
//			else if(lod > 3)
//			{
//				System.out.println("lod >3 "+geneName+",  " + group + "  has LOD " + lod);
//			}
//			
//			else if(lod > 2)
//			{
//				System.out.println("lod >2 "+geneName+", " + group + " has LOD " + lod);
//			}

			
		}
		

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
		// System.out.println(pval);

		
		return groupToPval;
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

	public HashMap<String, HashMap<String, Boolean>> getSampleLOF(LinkedHashMap<String, String> sequenceFeatureLocations) throws Exception
	{
		HashMap<String, HashMap<String, Boolean>> geneToSampleToLOF = new HashMap<String, HashMap<String, Boolean>>();

		Iterator<Entity> vcfIter = vcfReader.iterator();

		// keep track of the 'line' in which variant was seen
		int index = 0;
		int passedFilter = 0;

		vcfIterNext: while (vcfIter.hasNext())
		{
			index++;

			if (index % 1000 == 0)
			{
				System.out.println("Seen " + index + " variants..");
			}

			if (devMode && index % 10000 == 0)
			{
				System.out.println("DEV !! quitting..."); break;
			}

			Entity record = vcfIter.next();

			// System.out.println("looking at " + record.toString());
			
			String filter = record.getString("FILTER");
			if(!filter.equals("PASS"))
			{
				continue;
			}
			passedFilter++;

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

			multiAnn:
			for (int i = 0; i < multiAnn.length; i++)
			{
				String ann = multiAnn[i];
				String[] annSplit = ann.split("\\|", -1);
				String impact = annSplit[2];

				if (!impact.equals("HIGH"))
				{
					continue multiAnn; //FIXME: correct ??
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
						
						if(!sequenceFeatureLocations.containsKey(gene))
						{
							sequenceFeatureLocations.put(gene, chr + "\t" + pos);
						}

					}

				}

			}

		}
		
		System.out.println(passedFilter + " of " + index + " variants passed filter");

		//TODO: we don't have the genes that didn't match the criteria for at least 1 variant! e.g. for knockout we have ~500 genes...
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
	
	
	/**
	 * Create R script that plots image at given location using given dataset location
	 * 
	 * @param scriptLocation
	 * @param plotLocation
	 * @throws UnsupportedEncodingException
	 * @throws FileNotFoundException
	 */
	public void Rscript(String scriptLocation, String plotLocation, String dataframeLocation)
			throws FileNotFoundException, UnsupportedEncodingException
	{

		PrintWriter scriptPw = new PrintWriter(scriptLocation, "UTF-8");
		scriptPw.println("library(qqman)");
		scriptPw.println("hits <- as.data.frame(read.table(\"" + dataframeLocation
				+ "\", check.names=FALSE, header=TRUE, sep =\"\\t\", quote=\"\", as.is=TRUE))");
		for(String group : groups)
		{
			scriptPw.println("hits$P_"+group+" <- as.numeric(hits$P_"+group+")");
		}
		//preselect 1 group as "the P-value"
		scriptPw.println("hits$P <- as.numeric(hits$P_"+groups.toArray()[0].toString()+")");
		scriptPw.println("png(\"" + plotLocation + "\", res=200, width=1920, height=1080)");
		scriptPw.println("manhattan(hits, col = c(\"blue4\", \"orange3\"))");
		scriptPw.println("dev.off()");
		scriptPw.flush();
		scriptPw.close();

		RScriptExecutor r = new RScriptExecutor("/usr/bin/Rscript", null);
		ROutputHandler outputHandler = new StringROutputHandler();
		r.executeScript(new File(scriptLocation), outputHandler);

	}


}
