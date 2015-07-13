package org.molgenis.caddtlmapping;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Scanner;

import org.molgenis.data.AttributeMetaData;
import org.molgenis.data.Entity;
import org.molgenis.data.annotation.RepositoryAnnotator;
import org.molgenis.data.support.DefaultAttributeMetaData;
import org.molgenis.data.support.DefaultEntityMetaData;
import org.molgenis.data.vcf.VcfRepository;

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
	 * the index of candidates.
	 * 
	 * @param vcfFile
	 * @param exacFile
	 * @param mafThreshold
	 * @return
	 * @throws IOException
	 */
	public static HashMap<String, List<Integer>> getCandidateVariantsPerGene(File vcfFile, File exacFile,
			double mafThreshold) throws IOException
	{
		VcfRepository vcfRepo = new VcfRepository(vcfFile, "CADDTLMapping");
		Iterator<Entity> vcfRepoIter = vcfRepo.iterator();
		while (vcfRepoIter.hasNext())
		{
			Entity record = vcfRepoIter.next();
			System.out.println(record.toString());
		}
		return null;
	}
}
