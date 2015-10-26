package org.molgenis.calibratecadd;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

import org.molgenis.calibratecadd.misc.Step4_Helper;
import org.molgenis.calibratecadd.structs.EntityPlus;
import org.molgenis.calibratecadd.structs.ImpactRatios;
import org.molgenis.calibratecadd.structs.VariantIntersectResult;
import org.molgenis.data.Entity;
import org.molgenis.data.annotator.tabix.TabixVcfRepository;
import org.molgenis.data.vcf.VcfRepository;

public class Step4_MatchingVariantsFromExAC
{

	/**
	 * Uses:
	 * [0] file produced in step 3
	 * [1] ftp://ftp.broadinstitute.org/pub/ExAC_release/release0.3/ExAC.r0.3.sites.vep.vcf.gz (+ in the same folder ExAC.r0.3.sites.vep.vcf.gz.tbi )
	 * [2] output file
	 * 
	 * Example:
	 * E:\Data\clinvarcadd\clinvar.patho.fix.snpeff.vcf
	 * E:\Data\clinvarcadd\ExAC.r0.3.sites.vep.vcf.gz
	 * E:\Data\clinvarcadd\clinvar.patho.fix.snpeff.exac.vcf
	 * 
	 */
	public static void main(String[] args) throws Exception
	{
		Step4_MatchingVariantsFromExAC step4 = new Step4_MatchingVariantsFromExAC();
		step4.loadClinvarPatho(args[0]);
		step4.createMatchingExACsets(args[1]);
		step4.printVariantsToFile(args[2]);

	}

	HashMap<String, List<Entity>> clinvarPatho = new HashMap<String, List<Entity>>();
	HashMap<String, List<EntityPlus>> matchedExACvariants = new HashMap<String, List<EntityPlus>>();

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

	private void createMatchingExACsets(String exacLoc) throws Exception
	{
		System.out.println("loading matching exac variants..");

		int passedGenes = 0;
		int matchedVariants = 0;
		int droppedGenesClinVarTooFew = 0;
		int droppedGenesExACtooFew = 0;
		int droppedGenesNoMatchedVariants = 0;

		int index = 0;
		for (String gene : clinvarPatho.keySet())
		{
			/** DEV **/
//			if(!gene.equals("LOXHD1"))
//			{
//				continue;
//			}
			
			index++;
			if (index % 100 == 0)
			{
				System.out.println("Seen " + index + " of " + clinvarPatho.keySet().size() + " genes..");
			}

			String chrom = null;
			long leftMostPos = -1;
			long rightMostPos = -1;
			
			//we want at least 2 variants in order to get the interval for exac
			if(clinvarPatho.get(gene).size() < 2)
			{
				droppedGenesClinVarTooFew++;
				continue;
			}

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

			TabixVcfRepository tr = new TabixVcfRepository(new File(exacLoc), "exac");

			List<Entity> exacVariants = new ArrayList<Entity>();

			try
			{
				exacVariants = tr.query(chrom, leftMostPos, rightMostPos);
			}
			catch (java.lang.ArrayIndexOutOfBoundsException e)
			{
				// no chrom in tabix or so
			}

			System.out.println("\n#####\n");
			System.out.println(gene + " " + leftMostPos + " " + rightMostPos + " " + " has " + " " + exacVariants.size());

			if (exacVariants.size() > 0)
			{
				
				//found out: which variants are only in ExAC, which only in ClinVar, which in both
				VariantIntersectResult vir = Step4_Helper.intersectVariants(exacVariants, clinvarPatho.get(gene));
				
				System.out.println("VariantIntersectResult for '"+gene+"', clinvaronly: " + vir.inClinVarOnly.size() + ", exaconly: " + vir.inExACOnly.size() + ", both: " + vir.inBoth_exac.size());
				
				//calculate MAF for shared variants, and use them to filter the other ExAC variants
				//this way, we use the overlap to determine a fair cutoff for the 'assumed benign' variants
				//if we have nothing to go on, we will set this to 0 and only select singleton variants
				double medianMAF = Step4_Helper.calculateMedianMAF(vir.inBoth_exac);
				List<EntityPlus> exacFilteredByMAF = Step4_Helper.filterExACvariantsByMAF(vir.inExACOnly, medianMAF);
				
				System.out.println("exaconly filtered down to " + exacFilteredByMAF.size() + " variants using MAF " + medianMAF);

				//calculate impact ratios over all clinvar variants, and use them to 'shape' the remaining ExAC variants
				//they must become a set that looks just like the ClinVar variants, including same distribution of impact types
				ImpactRatios ir = Step4_Helper.calculateImpactRatios(vir.inClinVarOnly, vir.inBoth_clinvar);
				System.out.println("impact ratio from clinvar variants: " + ir.toString());
				List<EntityPlus> exacFilteredByMAFandImpact = Step4_Helper.shapeExACvariantsByImpactRatios(exacFilteredByMAF, ir);
				System.out.println("exaconly filtered down to " + exacFilteredByMAFandImpact.size() + " variants");
				
				if (exacFilteredByMAFandImpact.size() > 0)
				{
					passedGenes++;
					matchedExACvariants.put(gene, exacFilteredByMAFandImpact);
					matchedVariants += exacFilteredByMAFandImpact.size();
				}
				else
				{
					droppedGenesNoMatchedVariants++;
				}
				
				

			}
			else
			{
				droppedGenesExACtooFew++;
			}
			
			tr.close();
			
		}
		
		System.out.println();
		System.out.println("#### done ####");
		System.out.println();
		
		// oct 2015: 2638 pass, 960040 variants, 393 dropped
		System.out.println("passed genes (>0 properly matched interval exac variants): " + passedGenes);
		System.out.println("matched variants (total variants used for final calibration): " + matchedVariants);
		System.out.println("dropped genes (less than 2 clinvar variants): " + droppedGenesClinVarTooFew);
		System.out.println("dropped genes (2+ clinvar, but 0 interval exac variants): " + droppedGenesExACtooFew);
		System.out.println("dropped genes (2+ clinvar, but >0 interval exac variants, but 0 matched variants left after filtering): " + droppedGenesNoMatchedVariants);
		
	}
	
	
	private void printVariantsToFile(String file) throws FileNotFoundException
	{
		PrintWriter pw_variantInfo = new PrintWriter(file + ".info");
		PrintWriter pw_forCADD = new PrintWriter(file + ".cadd");
		
		pw_variantInfo.println( "gene" + "\t" + "chr" + "\t" + "pos" + "\t" + "ref" + "\t" + "alt" + "\t" + "group");
		
		for(String gene : clinvarPatho.keySet())
		{
			if(matchedExACvariants.containsKey(gene))
			{
				//print data from clinvarPatho and matchedExACvariants to file
				for(Entity variant : clinvarPatho.get(gene))
				{
					pw_forCADD.println(variant.getString("#CHROM") + "\t" + variant.getString("POS") + "\t" + "." + "\t" + variant.getString("REF") + "\t" + variant.getString("ALT"));
					pw_variantInfo.println(gene + "\t" + variant.getString("#CHROM") + "\t" + variant.getString("POS") + "\t" + variant.getString("REF") + "\t" + variant.getString("ALT") + "\t" + "PATHOGENIC");
				}
				for(EntityPlus variant : matchedExACvariants.get(gene))
				{
					pw_forCADD.println(variant.getE().getString("#CHROM") + "\t" + variant.getE().getString("POS") + "\t" + "." + "\t"+ variant.getE().getString("REF") + "\t" + variant.getKeyVal().get("ALT").toString());
					pw_variantInfo.println(gene + "\t" + variant.getE().getString("#CHROM") + "\t" + variant.getE().getString("POS") + "\t" + variant.getE().getString("REF") + "\t" + variant.getKeyVal().get("ALT").toString() + "\t" + "POPULATION");
				}
			}
		}
		
		
		
		
		pw_variantInfo.flush();
		pw_variantInfo.close();
		pw_forCADD.flush();
		pw_forCADD.close();
	}


}
