package org.molgenis.caddtlmapping.mannwhitney;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.UnsupportedEncodingException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Scanner;

import org.molgenis.data.Entity;
import org.molgenis.data.annotation.entity.impl.SnpEffAnnotator;
import org.molgenis.data.annotator.tabix.TabixReader;
import org.molgenis.data.annotator.tabix.TabixVcfRepository;
import org.molgenis.data.vcf.VcfRepository;
import org.molgenis.r.ROutputHandler;
import org.molgenis.r.RScriptExecutor;
import org.molgenis.r.StringROutputHandler;

import com.google.common.collect.ImmutableList;

public class Helper
{
	// private int secondPassIndex;
	// private int secondPassIndexForPopContrasting;
	private TabixVcfRepository patientTabixReader;
	private TabixReader exacTabixReader;
	private TabixReader caddTabixReader;

	public Helper(File vcfFile, File exacFile, File caddFile) throws IOException
	{
		// this.secondPassIndex = 0;
		// this.secondPassIndexForPopContrasting = 0;
		this.patientTabixReader = new TabixVcfRepository(vcfFile, "patientreader");
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
	 * Returns a map of 'sequence feature' to list of 'candidate variants' denoted by chr:pos-pos location in the Tabix
	 * file.
	 * 
	 * @param vcfFile
	 * @param exacFile
	 * @param mafThreshold
	 * @return
	 * @throws Exception
	 */
	public LinkedHashMap<String, List<String>> getCandidateVariantsPerGene(File vcfFile, double mafThreshold,
			Map<String, String> sequenceFeatureLocations) throws Exception
	{

		LinkedHashMap<String, List<String>> res = new LinkedHashMap<String, List<String>>();

		VcfRepository vcfRepo = new VcfRepository(vcfFile, "CADDTLMappingFirstPass");
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
			// if (index == 10000)
			// {
			// vcfRepo.close();
			// return res;
			// }
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
			if (filter == null || !filter.equals("PASS"))
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
				System.out.println("INDEVELOPMENT: for now, ignoring multi allelic input patient variants!! "
						+ record.toString());
				continue;
			}

			// filter by protein impact
			// TODO for now, hard-code only consider MODERATE and HIGH impact, but later on, make this optional or
			// configurable

			// multiple annotations? could even be multiple genes!!
			// example:
			// ANN=T|missense_variant|MODERATE|ANGPTL7|ANGPTL7|transcript|NM_021146.3|Coding|3/5|c.658C>T|p.Arg220Cys|949/2290|658/1041|220/346||,T|intron_variant|MODIFIER|MTOR|MTOR|transcript|NM_004958.3|Coding|28/57|c.4253+5498G>A||||||

			String[] multiAnn = record.getString("INFO_ANN").split(",");

			for (int i = 0; i < multiAnn.length; i++)
			{
				String ann = multiAnn[i];
				String[] annSplit = ann.split("\\|", -1);
				SnpEffAnnotator.Impact impact = Enum.valueOf(SnpEffAnnotator.Impact.class, annSplit[2]);
				if (!(impact.equals(SnpEffAnnotator.Impact.MODERATE) || impact.equals(SnpEffAnnotator.Impact.HIGH)))
				{
					continue;
				}

				// filter by MAF
				double exacMaf = 0.0;
				String exacVariant = getFromExAC(chr, pos, ref, alt);
				if (exacVariant != null)
				{
					String AFstring = getInfoFieldFromExACVariant(exacVariant, "AF", chr, pos, ref, alt);
					exacMaf = Double.parseDouble(AFstring);
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
					List<String> candLocsForGene = new ArrayList<String>();
					res.put(gene, candLocsForGene);
					// for location, simply store pos of the first candidate here
					sequenceFeatureLocations.put(gene, chr + "\t" + pos);
				}

				// store index of this candidate
				// add 0.1 for other alt allele, maybe even 0.2 for second alt allele, etc. TODO
				// btw NOT the same as multiple snpEff annotations
				res.get(gene).add(chr + ":" + pos + "-" + pos);
			}

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
	private String getInfoFieldFromExACVariant(String exacVariant, String infoFieldWanted, String chr, String pos,
			String ref, String alt) throws Exception
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
			if (infoField.startsWith(infoFieldWanted + "="))
			{
				String infoValues = infoField.replace(infoFieldWanted + "=", "");
				String[] infoValuesSplit = infoValues.split(",", -1);
				String infoValueString = infoValuesSplit[indexOfAlt];
				return infoValueString;
			}
		}

		throw new Exception("Could not extract " + infoFieldWanted + " from variant: " + exacVariant);
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
	private String getFromExAC(String chr, String pos, String ref, String alt) throws IOException
	{

		TabixReader.Iterator it = this.exacTabixReader.query(chr + ":" + pos + "-" + pos);
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

	// private String getCandidateVariant(String chrPosQuery) throws Exception
	// {
	// TabixReader.Iterator it = patientTabixReader.query(chrPosQuery);
	// //1:23000-24000
	// String askedChr = chrPosQuery.split(":",-1)[0];
	// String askedPos = chrPosQuery.split("-",-1)[1];
	// String next;
	// while (it != null && (next = it.next()) != null)
	// {
	// String[] nextSplit = next.split("\t", -1);
	// if(nextSplit[0].equals(askedChr) && nextSplit[1].equals(askedPos))
	// {
	// return next;
	// }
	// }
	// throw new Exception("Could not find " + chrPosQuery + " in patient vcf!");
	// }

	/**
	 * Get list of CADD scores for indices within patient VCF adjusted by inheritance mode
	 * 
	 * @param candList
	 * @param vcfFile
	 * @param caddFile
	 * @param inheritance
	 * @return
	 * @throws Exception
	 */
	public double[] getPatientCaddScores(List<String> candList, String inheritance,
			ArrayList<String> patientSampleIdList) throws Exception
	{

		// HACK 2 CHECK
		// patientSampleIdList = null;

		// Entity record;
		// int variantsProcessed = 0;
		ArrayList<Double> allCaddScores = new ArrayList<Double>();

		for (String candChrPos : candList)
		{
	//		ImmutableList<Entity> candVariantList = patientTabixReader.queryTabixSyntax(candChrPos); FIXME
			ImmutableList<Entity> candVariantList = null;
			if (candVariantList.size() > 1)
			{
				// TODO: throw new Exception ??
				// or iterate over multiple hits ??
				System.out.println("WARNING: More than 1 candidate variant in query results for " + candChrPos + "!!");
			}
			Entity candVariant = candVariantList.get(0);
			String chr = candVariant.getString("#CHROM");
			String pos = candVariant.getString("POS");
			String ref = candVariant.getString("REF");
			String alt = candVariant.getString("ALT");

			Double cadd = getCaddScore(chr, pos, ref, alt);

			// System.out.println("YES! cadd: " + cadd);

			if (cadd == null)
			{
				System.out.println("CADD null, have to skip this variant...");

			}

			else
			{

				// iterate over the genotypes
				Iterable<Entity> sampleEntities = candVariant.getEntities(VcfRepository.SAMPLES);
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
						// not a patient, ignore
					}
				}

			}

			// variantsProcessed++;
		}

		// System.out.println("list.size()=" + list.size() + ", variantsProcessed=" + variantsProcessed);
		// if (candList.size() == variantsProcessed)
		// {
		//
		// // System.out.println("returning: " + allCaddScores.toString());
		//
		// }

		return convertDoubles(allCaddScores);

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
			return res;
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
	 * @param candList
	 * @param vcfFile
	 * @param caddFile
	 * @param exacFile
	 * @param inheritance
	 * @return
	 * @throws Exception
	 */
	public double[] getPopulationCaddScores(List<String> candList, String inheritance) throws Exception
	{

		ArrayList<Double> allCaddScores = new ArrayList<Double>();

		for (String candChrPos : candList)
		{
		//	ImmutableList<Entity> candVariantList = patientTabixReader.queryTabixSyntax(candChrPos); FIXME
			ImmutableList<Entity> candVariantList = null;
			if (candVariantList.size() > 1)
			{
				// TODO: throw new Exception ??
				// or iterate over multiple hits ??
				System.out.println("WARNING: More than 1 candidate variant in query results for " + candChrPos + "!!");
			}
			Entity candVariant = candVariantList.get(0);
			String chr = candVariant.getString("#CHROM");
			String pos = candVariant.getString("POS");
			String ref = candVariant.getString("REF");
			String alt = candVariant.getString("ALT");

			Double cadd = getCaddScore(chr, pos, ref, alt);

			String exacVariant = getFromExAC(chr, pos, ref, alt);

			if (exacVariant != null && cadd != null)
			{
				List<Double> caddForGeno = getGenotypeCountsForExACVariantAndCombineWithCADD(exacVariant, inheritance,
						cadd, chr, pos, ref, alt);
				allCaddScores.addAll(caddForGeno);

			}
			else
			{
				// TODO
				// instead, we should look for a comparable variant nearby!! or try other CADD files!! etc
			}

		}

		return convertDoubles(allCaddScores);
	}

	private List<Double> getGenotypeCountsForExACVariantAndCombineWithCADD(String exacVariant, String inheritance,
			Double cadd, String chr, String pos, String ref, String alt) throws Exception
	{
		List<Double> caddForExacVariantGenotypes = new ArrayList<Double>();

		String AC_HetString = getInfoFieldFromExACVariant(exacVariant, "AC_Het", chr, pos, ref, alt);
		String AC_HomString = getInfoFieldFromExACVariant(exacVariant, "AC_Hom", chr, pos, ref, alt);

		// if (AC_Het == null || AC_Hom == null)
		// {
		// throw new Exception("AC_Het == null || AC_Hom == null for " + exacVariant);
		// }

		int AC_Het = Integer.parseInt(AC_HetString);
		int AC_Hom = Integer.parseInt(AC_HomString);

		// for hom alt, always add twice
		for (int i = 0; i < AC_Hom; i++)
		{
			caddForExacVariantGenotypes.add(cadd);
			caddForExacVariantGenotypes.add(cadd);
		}

		for (int i = 0; i < AC_Het; i++)
		{
			// for het, it depends on inheritance model
			if (inheritance.equals("additive"))
			{
				// add once for additive
				caddForExacVariantGenotypes.add(cadd);
			}
			else if (inheritance.equals("dominant"))
			{
				// add twice for dominant
				caddForExacVariantGenotypes.add(cadd);
				caddForExacVariantGenotypes.add(cadd);
			}
			else if (inheritance.equals("recessive"))
			{
				// dont add them for recessive!
			}
			else
			{
				throw new Exception("unknown inheritance!" + inheritance);
			}
		}

		return caddForExacVariantGenotypes;
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
		scriptPw.println("hits$P <- as.numeric(hits$P)");
		scriptPw.println("png(\"" + plotLocation + "\", res=200, width=1920, height=1080)");
		scriptPw.println("manhattan(hits, col = c(\"blue4\", \"orange3\"))");
		scriptPw.println("dev.off()");
		scriptPw.flush();
		scriptPw.close();

		RScriptExecutor r = new RScriptExecutor("/usr/bin/Rscript", null);
		ROutputHandler outputHandler = new StringROutputHandler();
		r.executeScript(new File(scriptLocation), outputHandler);

	}

	/**
	 * Check if string is a number
	 * 
	 * @param str
	 * @return
	 */
	private static boolean isNumeric(String str)
	{
		try
		{
			Double.parseDouble(str);
		}
		catch (NumberFormatException nfe)
		{
			return false;
		}
		return true;
	}

	/**
	 * Needed for the qqman R package.. X should be 23, Y becomes 24, MT becomes 25 though this is arbitrary
	 * 
	 * @param loc
	 * @return
	 */
	public static String changeChromLetterToNumber(String loc)
	{
		String[] locSplit = loc.split("\t", -1);
		String chr = locSplit[0];

		// if numeric, good, just return the input loc
		if (isNumeric(chr))
		{
			return loc;
		}

		String pos = locSplit[1];
		if (chr.equalsIgnoreCase("X"))
		{
			return "23" + "\t" + pos;
		}
		else if (chr.equalsIgnoreCase("Y"))
		{
			return "24" + "\t" + pos;
		}
		else if (chr.equalsIgnoreCase("MT"))
		{
			return "25" + "\t" + pos;
		}
		else
		{
			System.out.println("WARNING: could not change chromosome letter " + chr
					+ " to a number! R plot using qqman will fail");
			return loc;
		}
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

	/**
	 * Sort hashmap
	 * 
	 * @param passedMap
	 * @return
	 */
	public LinkedHashMap<String, Double> sortHashMapByValuesD(HashMap<String, Double> passedMap)
	{
		List<String> mapKeys = new ArrayList<String>(passedMap.keySet());
		List<Double> mapValues = new ArrayList<Double>(passedMap.values());
		Collections.sort(mapValues);
		Collections.sort(mapKeys);

		LinkedHashMap<String, Double> sortedMap = new LinkedHashMap<String, Double>();

		Iterator<Double> valueIt = mapValues.iterator();
		while (valueIt.hasNext())
		{
			Object val = valueIt.next();
			Iterator<String> keyIt = mapKeys.iterator();

			while (keyIt.hasNext())
			{
				Object key = keyIt.next();
				String comp1 = passedMap.get(key).toString();
				String comp2 = val.toString();

				if (comp1.equals(comp2))
				{
					passedMap.remove(key);
					mapKeys.remove(key);
					sortedMap.put((String) key, (Double) val);
					break;
				}

			}

		}
		return sortedMap;
	}
}
