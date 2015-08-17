package org.molgenis.caddtlmapping.binom;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Scanner;

import org.molgenis.data.annotator.tabix.TabixReader;
import org.molgenis.data.annotator.tabix.TabixVcfRepository;

public class Helper
{
	private TabixVcfRepository patientTabixReader;
	private TabixVcfRepository exacTabixReader;
	private List<TabixReader> caddTabixReader;

	public Helper(File vcfFile, File exacFile, File caddFolder) throws IOException
	{
		this.patientTabixReader = new TabixVcfRepository(vcfFile, "patients");
		System.out.println("Created TabixVcfRepository 'patients' on " + vcfFile.getName());
		
		this.exacTabixReader = new TabixVcfRepository(exacFile, "population");
		System.out.println("Created TabixVcfRepository 'population' on " + exacFile.getName());
		
		this.caddTabixReader = new ArrayList<TabixReader>();
		for(File f : caddFolder.listFiles())
		{
			if(f.getName().endsWith(".tsv.gz"))
			{
				caddTabixReader.add(new TabixReader(f.getAbsolutePath()));
				System.out.println("Created TabixReader on '"+f.getName()+"'");
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

	public HashMap<String, List<BinnedGenotypeCount>> getBinnedGenotypeCountsPerSequenceFeature(File exacFile, double[] bins,
			String inheritance) throws Exception
	{
		if(bins.length % 2 != 0)
		{
			throw new Exception("getBinnedGenotypeCountsPerSequenceFeature: bins.length % 2 != 0");
		}
		
		HashMap<String, List<BinnedGenotypeCount>> res = new HashMap<String, List<BinnedGenotypeCount>>();
		
		
		
		
		
		
		return res;
	}

	
}
