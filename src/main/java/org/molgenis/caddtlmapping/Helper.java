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
	 * @throws IOException
	 */
	public static HashMap<String, List<Double>> getCandidateVariantsPerGene(File vcfFile, File exacFile,
			double mafThreshold, HashMap<String, String> sequenceFeatureLocations) throws IOException
	{
		VcfRepository vcfRepo = new VcfRepository(vcfFile, "CADDTLMapping");
		Iterator<Entity> vcfRepoIter = vcfRepo.iterator();
		int index = 0;
		while (vcfRepoIter.hasNext())
		{
			if(index % 5000 == 0)
			{
				System.out.println("Seen " + index + " variants..");
			}
			
			Entity record = vcfRepoIter.next();
			// System.out.println(record.toString());
			index++;
		}
		vcfRepo.close();
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
