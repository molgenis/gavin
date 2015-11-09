package org.molgenis.calibratecadd;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import org.apache.commons.lang3.StringUtils;
import org.molgenis.calibratecadd.support.EntityPlus;
import org.molgenis.calibratecadd.support.ImpactRatios;
import org.molgenis.calibratecadd.support.Step4_Helper;
import org.molgenis.calibratecadd.support.VariantIntersectResult;
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

	public static String NA = "x";
	HashMap<String, List<Entity>> clinvarPatho = new HashMap<String, List<Entity>>();
	HashMap<String, List<EntityPlus>> matchedExACvariants = new HashMap<String, List<EntityPlus>>();
	HashMap<String, String> geneInfo = new HashMap<String, String>();

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
			index++;
			
			/** DEV **/
//			if(!gene.equals("MAGEL2"))
//			{
//				continue;
//			}
//			if(index % 100 == 0)
//			{
//				break;
//			}

			String chrom = null;
			long leftMostPos = -1;
			long rightMostPos = -1;
			
			//we want at least 2 variants in order to get the interval for exac
			//we expect clinvarPatho.get(gene).size() == 1, but lets print anyway to check
			if(clinvarPatho.get(gene).size() < 2)
			{
				droppedGenesClinVarTooFew++;
				geneInfo.put(gene, "N1" + "\t" + clinvarPatho.get(gene).get(0).getString("#CHROM") + "\t" + clinvarPatho.get(gene).get(0).getString("POS") + "\t" + clinvarPatho.get(gene).get(0).getString("POS") + "\t" + NA + "\t" + clinvarPatho.get(gene).size() + StringUtils.repeat("\t" + NA, 5) + "\t" + Step4_Helper.calculateImpactRatiosFromClinVar(clinvarPatho.get(gene)) + StringUtils.repeat("\t" + NA, 4));
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
			System.out.println(gene + " (" + index + " of " + clinvarPatho.keySet().size() + ") " + leftMostPos + " " + rightMostPos + " " + " has " + exacVariants.size());

			if (exacVariants.size() > 0)
			{
				
				//found out: which variants are only in ExAC, which only in ClinVar, which in both
				VariantIntersectResult vir = Step4_Helper.intersectVariants(exacVariants, clinvarPatho.get(gene));
				
				System.out.println("VariantIntersectResult for '"+gene+"', clinvaronly: " + vir.inClinVarOnly.size() + ", exaconly: " + vir.inExACOnly.size() + ", both: " + vir.inBoth_exac.size());
				
				//calculate MAF for shared variants, and use them to filter the other ExAC variants
				//this way, we use the overlap to determine a fair cutoff for the 'assumed benign' variants
				//if we have nothing to go on, we will set this to 0 and only select singleton variants
				double pathogenicMAF = Step4_Helper.calculatePathogenicMAF(vir.inBoth_exac, vir.inClinVarOnly.size());
				List<EntityPlus> exacFilteredByMAF = Step4_Helper.filterExACvariantsByMAF(vir.inExACOnly, pathogenicMAF);

				System.out.println("exaconly filtered down to " + exacFilteredByMAF.size() + " variants using pathogenic MAF " + pathogenicMAF);

				//calculate impact ratios over all clinvar variants, and use them to 'shape' the remaining ExAC variants
				//they must become a set that looks just like the ClinVar variants, including same distribution of impact types
				ImpactRatios pathoImpactRatio = Step4_Helper.calculateImpactRatios(Stream.concat(vir.inClinVarOnly.stream(), vir.inBoth_clinvar.stream()).collect(Collectors.toList()));
				String unfilteredExacImpactRatio = vir.inExACOnly.size() > 0 ? Step4_Helper.calculateImpactRatios(vir.inExACOnly).toString() : NA;
				
				if(exacFilteredByMAF.size() == 0)
				{
					geneInfo.put(gene, "T1" + "\t" + chrom + "\t" + leftMostPos + "\t" + rightMostPos + "\t" + vir.inExACOnly.size() + "\t" + clinvarPatho.get(gene).size() + "\t" + pathogenicMAF + "\t" + unfilteredExacImpactRatio + "\t" + pathoImpactRatio.toString() + StringUtils.repeat("\t" + NA, 4));
					continue;
				}
				
				List<EntityPlus> exacFilteredByMAFandImpact = Step4_Helper.shapeExACvariantsByImpactRatios(exacFilteredByMAF, pathoImpactRatio);
				System.out.println("exaconly filtered down to " + exacFilteredByMAFandImpact.size() + " variants");
				
				
				if (exacFilteredByMAFandImpact.size() > 0)
				{
					passedGenes++;
					matchedExACvariants.put(gene, exacFilteredByMAFandImpact);
					matchedVariants += exacFilteredByMAFandImpact.size();
					//impacts AFTER correction
					ImpactRatios MAFandImpactFilteredExacImpactRatio = Step4_Helper.calculateImpactRatios(exacFilteredByMAFandImpact);
					geneInfo.put(gene, "Cx" + "\t" + chrom + "\t" + leftMostPos + "\t" + rightMostPos + "\t" + exacFilteredByMAFandImpact.size() + "\t" + clinvarPatho.get(gene).size() + "\t" + pathogenicMAF + "\t" + unfilteredExacImpactRatio + "\t" + pathoImpactRatio.toString().toString() + "\t" + MAFandImpactFilteredExacImpactRatio.toString());
				}
				else
				{
					//impacts BEFORE correction (which set it to 0)
					ImpactRatios MAFfilteredExacImpactRatio = Step4_Helper.calculateImpactRatios(exacFilteredByMAF);
					String cat = Step4_Helper.determineImpactFilterCat(MAFfilteredExacImpactRatio, pathoImpactRatio, pathogenicMAF);
					geneInfo.put(gene, cat + "\t" + chrom + "\t" + leftMostPos + "\t" + rightMostPos + "\t" + exacFilteredByMAF.size() + "\t" + clinvarPatho.get(gene).size() + "\t" + pathogenicMAF + "\t" + unfilteredExacImpactRatio + "\t" + pathoImpactRatio.toString() + "\t" + MAFfilteredExacImpactRatio.toString());
					droppedGenesNoMatchedVariants++;
				}
				
				

			}
			else
			{
				droppedGenesExACtooFew++;
				//can happen: 2+ ClinVar variants on same position, so we query an interval of size 1... so be it
				geneInfo.put(gene, "N2" + "\t" + chrom + "\t" + leftMostPos + "\t" + rightMostPos + "\t" + 0 + "\t" + clinvarPatho.get(gene).size() + StringUtils.repeat("\t" + NA, 5)  + "\t" + Step4_Helper.calculateImpactRatiosFromClinVar(clinvarPatho.get(gene)) + StringUtils.repeat("\t" + NA, 4));
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
		PrintWriter pw_variantInfo = new PrintWriter(file + ".variants.tsv");
		PrintWriter pw_forCADD = new PrintWriter(file + ".cadd.tsv");
		PrintWriter pw_geneInfo = new PrintWriter(file + ".genes.tsv");
		
		pw_variantInfo.println( "gene" + "\t" + "chr" + "\t" + "pos" + "\t" + "ref" + "\t" + "alt" + "\t" + "group");
		pw_geneInfo.println( "Gene" + "\t" + "Category" + "\t" + "Chr" + "\t" + "Start" + "\t" + "End" + "\t" + "NrOfPopulationVariants" + "\t" + "NrOfPathogenicVariants" + "\t" + "PathoMAFThreshold" + "\t" + "PopImpactHighPerc" + "\t" + "PopImpactModeratePerc" + "\t" + "PopImpactLowPerc" + "\t" + "PopImpactModifierPerc" + "\t" + "PathoImpactHighPerc" + "\t" + "PathoImpactModeratePerc" + "\t" + "PathoImpactLowPerc" + "\t" + "PathoImpactModifierPerc" + "\t" + "PopImpactHighEq" + "\t" + "PopImpactModerateEq" + "\t" + "PopImpactLowEq" + "\t" + "PopImpactModifierEq");
		
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
			//replace "/" by "_" because R should not write output files with "/" in them, for obvious reasons.
			pw_geneInfo.println(gene.replace("/", "_") + "\t" + geneInfo.get(gene));
		}
		
		pw_variantInfo.flush();
		pw_variantInfo.close();
		
		pw_geneInfo.flush();
		pw_geneInfo.close();
		
		pw_forCADD.flush();
		pw_forCADD.close();
	}


}
