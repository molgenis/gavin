package org.molgenis.caddtlmapping;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
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
import org.molgenis.data.annotation.entity.impl.SnpEffAnnotator;
import org.molgenis.data.annotator.tabix.TabixReader;
import org.molgenis.data.vcf.VcfRepository;
import org.molgenis.r.ROutputHandler;
import org.molgenis.r.RScriptExecutor;
import org.molgenis.r.StringROutputHandler;

public class Helper
{
	private int secondPassIndex;
	private int secondPassIndexForPopContrasting;
	private TabixReader exacTabixReader;
	private TabixReader caddTabixReader;

	public Helper(File exacFile, File caddFile) throws IOException
	{
		this.secondPassIndex = 0;
		this.exacTabixReader = new TabixReader(exacFile.getAbsolutePath());
		this.caddTabixReader = new TabixReader(caddFile.getAbsolutePath());
	}

	/**
	 * Load sample id list
	 * 
	 * @param patientSampleIdList
	 * @return
	 * @throws FileNotFoundException
	 */
	public ArrayList<String> loadPatientSampleIdList(File patientSampleIdList) throws FileNotFoundException
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
	public HashMap<String, List<Double>> getCandidateVariantsPerGene(File vcfFile, double mafThreshold,
			HashMap<String, String> sequenceFeatureLocations) throws Exception
	{

		HashMap<String, List<Double>> res = new HashMap<String, List<Double>>();

		VcfRepository vcfRepo = new VcfRepository(vcfFile, "CADDTLMapping");
		Iterator<Entity> vcfRepoIter = vcfRepo.iterator();

		// keep track of the 'line' in which variant was seen
		int index = 0;

		while (vcfRepoIter.hasNext())
		{
			index++;

			if (index % 1000 == 0)
			{
				System.out.println("Seen " + index + " variants..");
			}

			// FIXME
			/** DEV **/
			if (index == 2000)
			{
				vcfRepo.close();
				return res;
			}
			/** DEV **/

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

			// TODO we also expect the VCF to be sorted !! need to check

			String chr = record.getString("#CHROM");
			String pos = record.getString("POS");
			String ref = record.getString("REF");
			String alt = record.getString("ALT");
			String[] altSplit = alt.split(",", -1);

			if (altSplit.length > 1)
			{
				// TODO
				// System.out.println("WARNING: for now, ignoring multi allelic input patient variants!! " +
				// record.toString());
				continue;
			}

			// filter by protein impact
			// TODO for now, hard-code only consider MODERATE and HIGH impact, but later on, make this optional or
			// configurable
			String[] annSplit = record.getString("INFO_ANN").split("\\|", -1);
			SnpEffAnnotator.Impact impact = Enum.valueOf(SnpEffAnnotator.Impact.class, annSplit[2]);
			if (!(impact.equals(SnpEffAnnotator.Impact.MODERATE) || impact.equals(SnpEffAnnotator.Impact.HIGH)))
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
			// add 0.1 for other alt, maybe even 0.2 for second alt, etc. TODO
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
	private double getMAFForExACVariant(String exacVariant, String chr, String pos, String ref, String alt)
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
	private String getFromExAC(TabixReader exacTabixReader, String chr, String pos, String ref, String alt)
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
	 * @throws Exception
	 */
	public double[] getPatientCaddScores(List<Double> list, Iterator<Entity> patientVcfIter, String inheritance,
			ArrayList<String> patientSampleIdList) throws Exception
	{

		// HACK 2 CHECK
		// patientSampleIdList = null;

		Entity record;
		int variantsProcessed = 0;
		ArrayList<Double> allCaddScores = new ArrayList<Double>();
		while (patientVcfIter.hasNext())
		{
			this.secondPassIndex++;

			record = patientVcfIter.next();

			if (list.contains(Math.floor(secondPassIndex)))
			{

				String chr = record.getString("#CHROM");
				String pos = record.getString("POS");
				String ref = record.getString("REF");
				String alt = record.getString("ALT");

				Double cadd = getCaddScore(chr, pos, ref, alt);

				System.out.println("YES! cadd: " + cadd);

				if (cadd == null)
				{
					System.out.println("CADD null, have to skip this variant...");

				}

				else
				{

					// iterate over the genotypes
					Iterable<Entity> sampleEntities = record.getEntities(VcfRepository.SAMPLES);
					for (Entity sample : sampleEntities)
					{

						if (patientSampleIdList == null)
						{
							// TODO in this case: consider EVERY genotype!
							// refactor.......

							List<Double> caddForGeno = combineGenotypeWithCadd(sample.get("GT").toString(), cadd,
									inheritance);
							allCaddScores.addAll(caddForGeno);

						}

						// TODO known ugly....
						// NAME='905957_T_100386'
						else if (patientSampleIdList.contains(sample.get("NAME").toString().split("_")[2]))
						{
							List<Double> caddForGeno = combineGenotypeWithCadd(sample.get("GT").toString(), cadd,
									inheritance);
							allCaddScores.addAll(caddForGeno);
						}

						else
						{

						}
					}

				}

				variantsProcessed++;
			}

			if (list.size() == variantsProcessed)
			{

				// System.out.println("returning: " + allCaddScores.toString());
				return convertDoubles(allCaddScores);
			}

		}

		throw new Exception(
				"getPatientCaddScores ran through entire patient VCF file without finalizing candidate variants, something went wrong");
	}

	/**
	 * Combine genotype with CADD score and inheritance model to output stuff
	 * 
	 * @param geno
	 * @param cadd
	 * @param inheritance
	 * @return
	 * @throws Exception
	 */
	private List<Double> combineGenotypeWithCadd(String genotype, double cadd, String inheritance) throws Exception
	{
		// System.out.println("combineGenotypeWithCadd " + genotype + " " + cadd + " " + inheritance);
		ArrayList<Double> res = new ArrayList<Double>(2);

		if (genotype.equals("./."))
		{
			// missing genotype... skip
			return null;
		}

		String[] genoSplit = genotype.split("/", -1);
		if (genoSplit.length != 2)
		{
			throw new Exception("Genotype not splittable: " + genotype);
		}

		if (!(genoSplit[0].equals("0") || genoSplit[0].equals("1")))
		{
			throw new Exception("allele" + genoSplit[0] + " is not 0 or 1, others currently not supported");
		}

		if (!(genoSplit[1].equals("0") || genoSplit[1].equals("1")))
		{
			throw new Exception("allele" + genoSplit[1] + " is not 0 or 1, others currently not supported");
		}

		int allele0 = Integer.parseInt(genoSplit[0]);
		int allele1 = Integer.parseInt(genoSplit[1]);

		if (allele0 + allele1 == 0)
		{
			return res;
		}

		else if (allele0 + allele1 == 2)
		{
			res.add(cadd);
			res.add(cadd);
			return res;
		}

		else if (allele0 + allele1 == 1)
		{
			if (inheritance.equals("dominant"))
			{
				res.add(cadd);
				res.add(cadd);
				return res;
			}
			else if (inheritance.equals("recessive"))
			{
				return res;
			}
			else if (inheritance.equals("additive"))
			{
				res.add(cadd);
				return res;
			}
			else
			{
				throw new Exception("Inheritance " + inheritance + " unknown");
			}
		}
		else
		{
			throw new Exception("allele0 + allele1 unknown value " + (allele0 + allele1));
		}

	}

	/**
	 * Get CADD score
	 * 
	 * @param chr
	 * @param pos
	 * @param ref
	 * @param alt
	 * @return
	 * @throws IOException
	 */
	private Double getCaddScore(String chr, String pos, String ref, String alt) throws IOException
	{
		TabixReader.Iterator it = caddTabixReader.query(chr + ":" + pos + "-" + pos);
		String next;
		while (it != null && (next = it.next()) != null)
		{
			String[] split = next.split("\t", -1);
			// 1 878331 C T 1.590426 13.80
			if (split[0].equals(chr) && split[1].equals(pos) && split[2].equals(ref) && split[3].equals(alt))
			{
				return Double.parseDouble(split[4]);
			}
		}
		return null;
	}

	/**
	 * Get a comparable set of variants from the population to contract against
	 * 
	 * @param list
	 * @param vcfFile
	 * @param caddFile
	 * @param exacFile
	 * @param inheritance
	 * @return
	 */
	public double[] getPopulationCaddScores(List<Double> list, Iterator<Entity> patientVcfIterForPopContrasting,
			String inheritance)
	{

		Entity record;
		int variantsProcessed = 0;
		ArrayList<Double> allCaddScores = new ArrayList<Double>();
		while (patientVcfIterForPopContrasting.hasNext())
		{
			this.secondPassIndexForPopContrasting++;

			record = patientVcfIterForPopContrasting.next();

			if (list.contains(Math.floor(secondPassIndexForPopContrasting)))
			{

				String chr = record.getString("#CHROM");
				String pos = record.getString("POS");
				String ref = record.getString("REF");
				String alt = record.getString("ALT");
				
				
				
			}
		}

		return null;
	}

	/**
	 * Helper function to sort LOD scores
	 * 
	 * @author jvelde
	 *
	 */
	public <K, V extends Comparable<? super V>> LinkedHashMap<String, Double> sortByValue(Map<String, Double> map)
	{
		LinkedHashMap<String, Double> result = new LinkedHashMap<>();
		Stream<Entry<String, Double>> st = map.entrySet().stream();
		st.sorted(Comparator.comparing(e -> e.getValue())).forEach(e -> result.put(e.getKey(), e.getValue()));
		return result;
	}

	private void Rscript()
	{
		/*
		 * png("~/Desktop/test.png"); plot(c(1,2,3,4,5)); dev.off()
		 */
		RScriptExecutor r = new RScriptExecutor("/usr/bin/Rscript");
		ROutputHandler outputHandler = new StringROutputHandler();
		r.executeScript(new File("/Users/jvelde/Desktop/test.R"), outputHandler);
	}

	/**
	 * Convert list to primitive array
	 * 
	 * @param doubles
	 * @return
	 */
	private double[] convertDoubles(List<Double> doubles)
	{
		double[] ret = new double[doubles.size()];
		for (int i = 0; i < ret.length; i++)
		{
			ret[i] = doubles.get(i).doubleValue();
		}
		return ret;
	}
}
