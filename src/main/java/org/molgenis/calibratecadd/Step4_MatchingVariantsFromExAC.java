package org.molgenis.calibratecadd;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

import org.molgenis.calibratecadd.structs.ClinVarVariant;
import org.molgenis.data.Entity;
import org.molgenis.data.annotator.tabix.TabixReader;
import org.molgenis.data.annotator.tabix.TabixReader.Iterator;
import org.molgenis.data.vcf.VcfRepository;

public class Step4_MatchingVariantsFromExAC
{

	/**
	 * Uses:
	 * [0] file produced in step 3
	 * [1] ftp://ftp.broadinstitute.org/pub/ExAC_release/release0.3/ExAC.r0.3.sites.vep.vcf.gz (+ in the same folder ExAC.r0.3.sites.vep.vcf.gz.tbi )
	 * [2] output file
	 */
	public static void main(String[] args) throws Exception
	{
		Step4_MatchingVariantsFromExAC step4 = new Step4_MatchingVariantsFromExAC();
		step4.loadClinvarPatho(args[0]);
		step4.getGeneIntervalVariantsFromExAC(args[1]);

	}

	HashMap<String, List<Entity>> clinvarPatho = new HashMap<String, List<Entity>>();
	
	HashMap<String, List<String>> exacIntervals = new HashMap<String, List<String>>();

	private void loadClinvarPatho(String clinvarPathoLoc) throws Exception
	{
		System.out.println("loading clinvar pathogenic variants from " + clinvarPathoLoc + " ..");
		VcfRepository vcfRepo = new VcfRepository(new File(clinvarPathoLoc), "clinvar");

		java.util.Iterator<Entity> vcfRepoIter = vcfRepo.iterator();

		while (vcfRepoIter.hasNext())
		{
			Entity record = vcfRepoIter.next();
			if (record.getString("ANN") == null)
			{
				vcfRepo.close();
				throw new Exception("Please annotate the VCF with a recent snpEff version! ANN field not found for "
						+ record.toString());
			}

			if (record.getString(Step1_GetClinVarPathogenic.CLINVAR_INFO) == null)
			{
				vcfRepo.close();
				throw new Exception("Did you create this VCF running Step1? " + Step1_GetClinVarPathogenic.CLINVAR_INFO
						+ " field not found for " + record.toString());
			}

			String geneAccordingToClinVar = record.getString(Step1_GetClinVarPathogenic.CLINVAR_INFO).split("\\|", -1)[1];
			String geneAccordingToSnpEff = record.getString("ANN").split("\\|")[3];

			// if snpeff has no gene, use the clinvar one (difference solved!)
			if (geneAccordingToSnpEff.isEmpty())
			{
				geneAccordingToSnpEff = geneAccordingToClinVar;
			}
			// if snpeff gene contains the clinvar one, use the clinvar one (difference solved!)
			else if (geneAccordingToSnpEff.contains(geneAccordingToClinVar))
			{
				geneAccordingToSnpEff = geneAccordingToClinVar;
			}
			// if clinvar gene contains the snpeff one, use the snpeff one (difference solved!)
			else if (geneAccordingToSnpEff.contains(geneAccordingToClinVar))
			{
				geneAccordingToClinVar = geneAccordingToSnpEff;
			}

			// if there is still a difference, use both.. we have about 120 of these (oct 2015)

			if (!geneAccordingToClinVar.equals(geneAccordingToSnpEff))
			{
				// System.out.println("Discrepancy detected! geneAccordingToClinVar="+geneAccordingToClinVar + ", " +
				// "geneAccordingToSnpEff=" +geneAccordingToSnpEff);
				geneAccordingToClinVar = geneAccordingToClinVar + "/" + geneAccordingToSnpEff;
			}

			if (clinvarPatho.containsKey(geneAccordingToClinVar))
			{
				clinvarPatho.get(geneAccordingToClinVar).add(record);
			}
			else
			{
				List<Entity> variants = new ArrayList<Entity>();
				variants.add(record);
				clinvarPatho.put(geneAccordingToClinVar, variants);
			}
		}
		vcfRepo.close();
	}

	private void getGeneIntervalVariantsFromExAC(String exacLoc) throws IOException
	{
		System.out.println("loading matching exac variants..");

		int passedGenes = 0;
		int passedVariants = 0;
		int matchedVariants = 0;
		int droppedGenes = 0;
		
		int index = 0;
		for (String gene : clinvarPatho.keySet())
		{
			index++;
			if (index % 100 == 0)
			{
				System.out.println("Seen " + index + " of "+clinvarPatho.keySet().size()+" genes..");
			}
			
			
			String chrom = null;
			long leftMostPos = -1;
			long rightMostPos = -1;

			for (Entity cvv : clinvarPatho.get(gene))
			{
				long pos = cvv.getLong("POS");
				String cvvchrom = cvv.getString("#CHROM");
				if (pos > rightMostPos)
				{
					rightMostPos = pos;
				}

				if (pos < leftMostPos || leftMostPos == -1)
				{
					leftMostPos = pos;
				}
				if (chrom == null)
				{
					chrom = cvvchrom;
				}
			}

			TabixReader tr = new TabixReader(exacLoc);
			Iterator it = null;
		
			try{
				it = tr.query(chrom + ":" + leftMostPos + "-" + rightMostPos);
			}
			catch(java.lang.ArrayIndexOutOfBoundsException e)
			{
				//chromosome not found, e.g. MT
			}
			
			int count = 0;
			String next = null;
			ArrayList<String> exacLines = new ArrayList<String>();
			while (it != null && (next = it.next()) != null)
			{
				// System.out.println(next);
				exacLines.add(next);
				count++;
			}

	//		System.out.println(gene + " " + leftMostPos + " " + rightMostPos + " " + " has " + " " + count);
			
			if(count > 0)
			{
				passedGenes++;
				passedVariants += exacLines.size();
				
			//	java.lang.OutOfMemoryError: Java heap space
			//	exacIntervals.put(gene, exacLines);
				
				List<String> matchingExAC = createMatchingExACset(clinvarPatho.get(gene), exacLines);
				
				matchedVariants += matchingExAC.size();
				
			}
			else
			{
				droppedGenes++;
			}

		}
		
		//oct 2015: 2638 pass, 960040 variants, 393 dropped
		System.out.println("passed genes (>0 interval exac variants): " + passedGenes);
		System.out.println("passed variants (total count in passed genes): " + passedVariants);
		System.out.println("dropped genes (0 interval exac variants): " + droppedGenes);

	}
	
	
	
	private List<String> createMatchingExACset(List<Entity> clinvar, List<String> exac)
	{
		int nrHighImpactClinVar = 0;
		
		List<String> res = new ArrayList<String>();
		
		return res;
	}

}
