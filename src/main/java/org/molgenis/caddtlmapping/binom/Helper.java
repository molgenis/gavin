package org.molgenis.caddtlmapping.binom;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Scanner;
import java.util.Set;

import org.molgenis.caddtlmapping.binom.structs.Bin;
import org.molgenis.caddtlmapping.binom.structs.BinnedGenotypeCount;
import org.molgenis.data.Entity;
import org.molgenis.data.annotator.tabix.TabixReader;
import org.molgenis.data.annotator.tabix.TabixVcfRepository;
import org.molgenis.data.vcf.VcfRepository;

public class Helper
{
	private TabixVcfRepository patientTabixReader;
	private TabixVcfRepository exacTabixReader;
	private List<TabixReader> caddTabixReader;

	private boolean dev = false; // FIXME remove for final version

	public Helper(File vcfFile, File exacFile, File caddFolder) throws IOException
	{
		this.patientTabixReader = new TabixVcfRepository(vcfFile, "patients");
		System.out.println("Created TabixVcfRepository 'patients' on " + vcfFile.getName());

		this.exacTabixReader = new TabixVcfRepository(exacFile, "population");
		System.out.println("Created TabixVcfRepository 'population' on " + exacFile.getName());

		this.caddTabixReader = new ArrayList<TabixReader>();
		for (File f : caddFolder.listFiles())
		{
			if (f.getName().endsWith(".tsv.gz"))
			{
				caddTabixReader.add(new TabixReader(f.getAbsolutePath()));
				System.out.println("Created TabixReader on '" + f.getName() + "'");
			}

		}
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

	public HashMap<String, List<BinnedGenotypeCount>> getBinnedPatientGenotypeCountsPerSequenceFeature(Bin[] bins,
			String inheritance, ArrayList<String> patientSampleIdList) throws Exception
	{
		HashMap<String, List<BinnedGenotypeCount>> res = new HashMap<String, List<BinnedGenotypeCount>>();

		Iterator<Entity> patientRepoIter = patientTabixReader.iterator();

		// keep track of the 'line' in which variant was seen
		int index = 0;

		while (patientRepoIter.hasNext())
		{
			index++;

			if (index % 10000 == 0)
			{
				System.out.println("Seen " + index + " patient variants..");
			}

			Entity record = patientRepoIter.next();

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

			Set<String> partOfFeature = new HashSet<String>();

			String[] multiAnn = record.getString("INFO_ANN").split(",");

			for (int i = 0; i < multiAnn.length; i++)
			{
				String ann = multiAnn[i];
				String[] annSplit = ann.split("\\|", -1);
				String gene = annSplit[3];
				if (gene.isEmpty())
				{
					throw new Exception("reminder: gene can be empty :) check for it");
				}
				partOfFeature.add(gene);
			}

			Double cadd = null;
			try
			{
				cadd = getCaddScore(chr, pos, ref, alt);
			}
			catch (java.lang.ArrayIndexOutOfBoundsException e)
			{
				System.out.println("caught  java.lang.ArrayIndexOutOfBoundsException for " + record.toString());
			}
			if (cadd == null)
			{
				// TODO
				// System.out.println("No CADD score found! skipping...");
				continue;
			}

			int totalGenotypes = 0;
			int actingGenotypes = 0;

			// iterate over the genotypes
			Iterable<Entity> sampleEntities = record.getEntities(VcfRepository.SAMPLES);
			for (Entity sample : sampleEntities)
			{

				String genotype = sample.get("GT").toString();

				if (genotype.equals("./."))
				{
					continue;
				}

				totalGenotypes++;

				if (patientSampleIdList == null)
				{
					// TODO in this case: consider EVERY genotype!
					// refactor.......

					actingGenotypes += getActingGenotype(genotype, inheritance);

				}

				// TODO known ugly....
				// NAME='905957_T_100386'
				else if (patientSampleIdList.contains(sample.get("NAME").toString().split("_")[2]))
				{

					actingGenotypes += getActingGenotype(genotype, inheritance);

				}

				else
				{
					// not a patient, ignore
				}
			}

			addCountsToCaddBins(partOfFeature, bins, cadd, res, actingGenotypes, totalGenotypes);

		}

		return res;
	}

	// TODO only works for 0 or 1 based genotypes!! does not support 1/2 etc !!
	public double getActingGenotype(String genotype, String inheritance) throws Exception
	{

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
			return 0;
		}

		else if (allele0 + allele1 == 2)
		{
			return 1;
		}

		else if (allele0 + allele1 == 1)
		{
			if (inheritance.equals("dominant"))
			{
				return 1;
			}
			else if (inheritance.equals("recessive"))
			{
				return 0;
			}
			else if (inheritance.equals("additive"))
			{
				return 0.5;
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

	public HashMap<String, List<BinnedGenotypeCount>> getBinnedGenotypeCountsPerSequenceFeature(Bin[] bins,
			String inheritance) throws Exception
	{

		HashMap<String, List<BinnedGenotypeCount>> res = new HashMap<String, List<BinnedGenotypeCount>>();

		Iterator<Entity> exacRepoIter = exacTabixReader.iterator();

		// keep track of the 'line' in which variant was seen
		int index = 0;

		while (exacRepoIter.hasNext())
		{
			index++;

			if (index % 10000 == 0)
			{
				System.out.println("Seen " + index + " ExAC variants..");
			}

			if (dev)
			{
				/** DEV **/
				if (index == 100000)
				{
					exacTabixReader.close();
					return res;
				}

			}

			Entity record = exacRepoIter.next();

			String chr = record.getString("#CHROM");
			String pos = record.getString("POS");
			String ref = record.getString("REF");
			String alt = record.getString("ALT");

			String[] altSplit = alt.split(",", -1);

			if (altSplit.length > 1)
			{
				// TODO
				// System.out.println("WARNING: skipping multi allelic for now..");
				continue;
			}

			Double cadd = null;
			try
			{
				cadd = getCaddScore(chr, pos, ref, alt);
			}
			catch (ArrayIndexOutOfBoundsException e)
			{
				System.out.println("caught java.lang.ArrayIndexOutOfBoundsException for " + record.toString() + " !!");
				continue;
			}

			if (cadd == null)
			{
				// TODO
				// System.out.println("No CADD score found! skipping...");
				continue;
			}

			// System.out.println(record.toString());

			Set<String> partOfFeature = new HashSet<String>();
			if (record.getString("INFO_CSQ") != null)
			{
				String[] info_csq_split = record.getString("INFO_CSQ").toString().split(",", -1);
				for (String info_csq : info_csq_split)
				{
					String sequenceFeature = info_csq.split("\\|", -1)[14];

					if (!sequenceFeature.isEmpty())
					{
						partOfFeature.add(sequenceFeature);
					}
					// System.out.println(sequenceFeature);
				}
			}

			// TODO: X/Y hemi ? calculate nr of genotypes etc
			record.getString("INFO_AC_Hemi");

			String AC_HetString = record.getString("INFO_AC_Het");
			String AC_HomString = record.getString("INFO_AC_Hom");
			String AN_AdjString = record.getString("INFO_AN_Adj");

			int AC_Het = Integer.parseInt(AC_HetString);
			int AC_Hom = Integer.parseInt(AC_HomString);
			int AN_Adj = Integer.parseInt(AN_AdjString);

			double totalGenotypes = AN_Adj / 2.0;
			double actingGenotypes = 0;
			if (inheritance.equals("dominant"))
			{
				actingGenotypes = AC_Het + AC_Hom;
			}
			else if (inheritance.equals("additive"))
			{
				actingGenotypes = (AC_Het / 2.0) + AC_Hom; // FIXME: 11/2 = 5, so should be ok ?
			}
			else if (inheritance.equals("recessive"))
			{
				actingGenotypes = AC_Hom;
			}
			else
			{
				throw new Exception("unknown mode of inheritance " + inheritance);
			}

			addCountsToCaddBins(partOfFeature, bins, cadd, res, actingGenotypes, totalGenotypes);

		}

		return res;
	}

	private void addCountsToCaddBins(Set<String> partOfFeature, Bin[] bins, double cadd,
			HashMap<String, List<BinnedGenotypeCount>> res, double actingGenotypes, double totalGenotypes)
	{
		for (String sequenceFeature : partOfFeature)
		{
			for (Bin bin : bins)
			{
				if (cadd > bin.lower && cadd < bin.upper)
				{

					if (res.containsKey(sequenceFeature))
					{
						boolean updateExistingBGC = false;

						for (BinnedGenotypeCount existingBGC : res.get(sequenceFeature))
						{
							if (existingBGC.bin == bin)
							{
								// System.out.println("BINS ARE EQUAL");
								existingBGC.actingGenotypes += actingGenotypes;
								existingBGC.totalGenotypes += totalGenotypes;
								updateExistingBGC = true;
								// System.out.println("for existing '" + sequenceFeature +
								// "' and existing bin, updated to " + existingBGC);
							}
						}

						if (!updateExistingBGC)
						{
							BinnedGenotypeCount bgc = new BinnedGenotypeCount(bin, (int) actingGenotypes,
									(int) totalGenotypes);
							// System.out.println("for existing '" + sequenceFeature + "', added " + bgc);
							res.get(sequenceFeature).add(bgc);
						}

					}
					else
					{
						ArrayList<BinnedGenotypeCount> bgc = new ArrayList<BinnedGenotypeCount>();
						bgc.add(new BinnedGenotypeCount(bin, (int) actingGenotypes, (int) totalGenotypes));
						res.put(sequenceFeature, bgc);
						// System.out.println("for new '" + sequenceFeature + "' created " + bgc);
					}

				}
			}
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
		for (TabixReader t : caddTabixReader)
		{
			TabixReader.Iterator it = t.query(chr + ":" + pos + "-" + pos);
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

		}
		return null;
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
