package org.molgenis.caddtlmapping;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Scanner;
import java.util.stream.Stream;

import org.molgenis.data.Entity;
import org.molgenis.data.annotation.entity.impl.SnpEffServiceAnnotator;
import org.molgenis.data.annotator.tabix.TabixReader;
import org.molgenis.data.vcf.VcfRepository;
import org.molgenis.r.ROutputHandler;
import org.molgenis.r.RScriptExecutor;
import org.molgenis.r.StringROutputHandler;

public class Helper
{

	/**
	 * Load sample id list
	 * 
	 * @param patientSampleIdList
	 * @return
	 * @throws FileNotFoundException
	 */
	public static ArrayList<String> loadPatientSampleIdList(File patientSampleIdList) throws FileNotFoundException
	{
		// load the sample identifiers of the patients within the VCF file
		// (we are not interested non-disease samples such as parents of trios etc)
		Scanner s = new Scanner(patientSampleIdList);
		ArrayList<String> sampleIds = new ArrayList<String>();
		while (s.hasNextLine())
		{
			sampleIds.add(s.nextLine());
		}
		s.close();
		return sampleIds;
	}

	/**
	 * Returns a map of sequence feature to list of candidate variant 'line numbers'. We simply count along and store
	 * the index of candidates. We use doubles to deal with multi allelic variants. For example, usually an index could
	 * be 154.0, but if there are multiple variants, they are encoded as 154.1 for the first alternative allele, and so
	 * on. Also store the (approximate) locations of the sequence features in a map.
	 * 
	 * @param vcfFile
	 * @param exacFile
	 * @param mafThreshold
	 * @return
	 * @throws Exception
	 */
	public static HashMap<String, List<Double>> getCandidateVariantsPerGene(File vcfFile, File exacFile,
			double mafThreshold, HashMap<String, String> sequenceFeatureLocations) throws Exception
	{
		TabixReader exacTabixReader = new TabixReader(exacFile.getAbsolutePath());

		HashMap<String, List<Double>> res = new HashMap<String, List<Double>>();

		VcfRepository vcfRepo = new VcfRepository(vcfFile, "CADDTLMapping");
		Iterator<Entity> vcfRepoIter = vcfRepo.iterator();

		// keep track of the 'line' in which variant was seen
		int index = 0;

		while (vcfRepoIter.hasNext())
		{
			if (index % 1000 == 0)
			{
				System.out.println("Seen " + index + " variants..");
			}
			index++;

			Entity record = vcfRepoIter.next();

			if (record.getString("INFO_ANN") == null)
			{
				vcfRepo.close();
				throw new Exception("Please annotate the VCF with a recent snpEff version! ANN field not found for "
						+ record.toString());
			}

			// filter by quality
			String filter = record.getString("FILTER");
			if (!filter.equals("PASS"))
			{
				continue;
			}

			String chr = record.getString("#CHROM");
			String pos = record.getString("POS");
			String ref = record.getString("REF");
			String alt = record.getString("ALT");
			String[] altSplit = alt.split(",", -1);

			if (altSplit.length > 1)
			{
				// TODO
				// System.out.println("WARNING: for now, ignoring multi allelic input patient variants!! " + record.toString());
				continue;
			}

			// filter by protein impact
			// TODO for now, hard-code only consider MODERATE and HIGH impact, but later on, make this optional or
			// configurable
			String[] annSplit = record.getString("INFO_ANN").split("\\|", -1);
			SnpEffServiceAnnotator.Impact impact = Enum.valueOf(SnpEffServiceAnnotator.Impact.class, annSplit[2]);
			if (!(impact.equals(SnpEffServiceAnnotator.Impact.MODERATE) || impact
					.equals(SnpEffServiceAnnotator.Impact.HIGH)))
			{
				continue;
			}

			// filter by MAF
			double exacMaf = 0.0;
			String exacVariant = getFromExAC(exacTabixReader, chr, pos, ref, alt);
			if (exacVariant != null)
			{
				exacMaf = getMAFForExACVariant(exacVariant, chr, pos, ref, alt);
			}

			// filter by MAF
			if (exacMaf > mafThreshold)
			{
				continue;
			}

			// all passed!
			String gene = annSplit[3];

			if (!res.containsKey(gene))
			{
				List<Double> indicesForGene = new ArrayList<Double>();
				res.put(gene, indicesForGene);
				// for location, simply store pos of the first candidate here
				sequenceFeatureLocations.put(gene, pos);
			}

			// store index of this candidate
			res.get(gene).add(index + 0.0);
			
		}
		vcfRepo.close();
		return res;
	}

	/**
	 * Parse variant from ExAC to get the MAF out eg. AF=1.209e-03; Must deal with multi allelic variants Crazy example:
	 * 15 66641732 has AF=0.112,2.471e-05,1.648e-05;
	 * 
	 * @param exacVariant
	 * @return
	 * @throws Exception
	 */
	private static double getMAFForExACVariant(String exacVariant, String chr, String pos, String ref, String alt)
			throws Exception
	{
		String[] split = exacVariant.split("\t", -1);

		// double check on chr, pos and ref
		if (!split[0].equals(chr) || !split[1].equals(pos) || !split[3].equals(ref))
		{
			throw new Exception("Chr=" + chr + ", pos=" + pos + " or ref=" + ref + " did NOT match the input variant: "
					+ exacVariant);
		}

		int indexOfAlt = 0;
		String[] altSplit = split[4].split(",", -1);
		boolean altMatch = false;
		for (int i = 0; i < altSplit.length; i++)
		{
			String altInExAc = altSplit[i];
			if (alt.equals(altInExAc))
			{
				indexOfAlt = i;
				altMatch = true;
			}
		}

		if (!altMatch)
		{
			throw new Exception("Alt=" + alt + " did not match input variant: " + exacVariant);
		}

		String[] infoSplit = split[7].split(";", -1);
		for (String infoField : infoSplit)
		{
			if (infoField.startsWith("AF="))
			{
				String mafValues = infoField.replace("AF=", "");
				String[] mafValuesSplit = mafValues.split(",", -1);
				String mafString = mafValuesSplit[indexOfAlt];
				double maf = Double.parseDouble(mafString);
				return maf;
			}
		}

		throw new Exception("Could not extract MAF from variant: " + exacVariant);
	}

	/**
	 * Get variant from ExAC, null if not found Check on exact match of chrom, pos, ref, and alt If multiple alts, then
	 * the list must include the asked for alt May return a line with multiple alt alleles
	 * 
	 * @param chrom
	 * @param pos
	 * @param ref
	 * @param alt
	 * @return
	 * @throws IOException
	 */
	private static String getFromExAC(TabixReader exacTabixReader, String chr, String pos, String ref, String alt)
			throws IOException
	{

		TabixReader.Iterator it = exacTabixReader.query(chr + ":" + pos + "-" + pos);
		String next;
		while (it != null && (next = it.next()) != null)
		{
			String[] split = next.split("\t", -1);

			// tabix may be fuzzy on the true positions... check it first
			if (!split[1].equals(pos))
			{
				continue;
			}

			// check ref
			if (!split[3].equals(ref))
			{
				// TODO seem to be a lot of these..
				// System.out.println("WARNING: found request variant in ExAC at chr=" + chr + ", pos=" + pos
				// + ", but ref did not match! (" + ref + " wanted, " + split[3] + " returned)");
				continue;
			}

			// check alt
			String[] altSplit = split[4].split(",", -1);
			boolean altMatch = false;
			for (String altInExAc : altSplit)
			{
				if (alt.equals(altInExAc))
				{
					altMatch = true;
				}
			}

			if (!altMatch)
			{
				// TODO seem to be a lot of these..
				// System.out.println("WARNING: found request variant in ExAC at chr=" + chr + ", pos=" + pos + ", ref="
				// + ref + ", but alt did not match! (" + alt + " wanted, possible alts returned: " + split[4] + ")");
				continue;
			}

			return next;
		}

		return null;
	}

	/**
	 * Get list of CADD scores for indices within patient VCF adjusted by inheritance mode
	 * 
	 * @param list
	 * @param vcfFile
	 * @param caddFile
	 * @param inheritance
	 * @return
	 */
	public static double[] getPatientCaddScores(List<Double> list, File vcfFile, File caddFile, String inheritance,
			ArrayList<String> patientSampleIdList)
	{
		// TODO Auto-generated method stub
		return null;
	}

	/**
	 * 
	 * @param list
	 * @param vcfFile
	 * @param caddFile
	 * @param exacFile
	 * @param inheritance
	 * @return
	 */
	public static double[] getPopulationCaddScores(List<Double> list, File vcfFile, File caddFile, File exacFile,
			String inheritance)
	{
		// TODO Auto-generated method stub
		return null;
	}

	/**
	 * Helper function to sort LOD scores
	 * 
	 * @author jvelde
	 *
	 */
	public static <K, V extends Comparable<? super V>> LinkedHashMap<String, Double> sortByValue(Map<String, Double> map)
	{
		LinkedHashMap<String, Double> result = new LinkedHashMap<>();
		Stream<Entry<String, Double>> st = map.entrySet().stream();
		st.sorted(Comparator.comparing(e -> e.getValue())).forEach(e -> result.put(e.getKey(), e.getValue()));
		return result;
	}

	public void Rscript()
	{
		/*
		 * png("~/Desktop/test.png"); plot(c(1,2,3,4,5)); dev.off()
		 */
		RScriptExecutor r = new RScriptExecutor("/usr/bin/Rscript");
		ROutputHandler outputHandler = new StringROutputHandler();
		r.executeScript(new File("/Users/jvelde/Desktop/test.R"), outputHandler);
	}
}
