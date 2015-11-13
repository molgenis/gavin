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

	public static String NA = "";
	HashMap<String, List<Entity>> clinvarPatho = new HashMap<String, List<Entity>>();
	HashMap<String, List<EntityPlus>> matchedExACvariants = new HashMap<String, List<EntityPlus>>();
	HashMap<String, String> geneInfo = new HashMap<String, String>();
	
	// keep track of which 'ANN' field index was used to match the gene to ClinVar symbol
	//e.g. if the matched gene is the second ANN field, we need to use that impact in the analysis
	HashMap<String, Integer> variantToNonZeroSnpEffGeneIndex = new HashMap<String, Integer>();

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
			
			//sanity check
			if(geneAccordingToClinVar.length() == 0)
			{
				vcfRepo.close();
				throw new Exception("geneAccordingToClinVar length 0");
			}
			
			//add all gene symbols provided by SnpEff to a list
			List<String> genesAccordingToSnpEff = new ArrayList<String>();
			String[] annSplit = record.getString("ANN").split(",", -1);
			for(String ann : annSplit)
			{
				String symbol = ann.split("\\|")[3];
				if(!symbol.isEmpty())
				{
					genesAccordingToSnpEff.add(symbol);
				}
			}
			
			String finalGeneSymbol = null;

			// if snpeff has no genes, use the clinvar one (difference solved!)
			// happens a lot for mitochondrial genes (e.g. MT-ND6, MT-CYB, MT-RNR1 etc)
			if (genesAccordingToSnpEff.isEmpty())
			{
				finalGeneSymbol = geneAccordingToClinVar;
			}

			String chrPosRefAlt = record.getString("#CHROM") + "_" + record.getString("POS") + "_" + record.getString("REF") + "_" + record.getString("ALT");
			
			// no joy? then we check all SnpEff symbols vs ClinVar
			if(finalGeneSymbol == null)
			{
				int index = 0;
				for(String snpEffGene : genesAccordingToSnpEff)
				{
					// if this snpeff gene equals the clinvar one, we're done
					if (snpEffGene.equals(geneAccordingToClinVar))
					{
						finalGeneSymbol = geneAccordingToClinVar;
					}

					// sometimes, a SnpEff symbol 'contains' the ClinVar symbol
					// happens quite often, e.g. TTN-AS1 contains TTN, INS-IGF2 contains INS etc.
					// if this happens, assign the ClinVar symbol (ie. the smaller one) and we're done
					else if (snpEffGene.contains(geneAccordingToClinVar))
					{
						finalGeneSymbol = geneAccordingToClinVar;
					}
					
					// opposite scenario: clinvar gene contains the snpeff one, use the snpeff one
					// seems to happen once: ClinVar CORO7-PAM16 contains SnpEff gene PAM16
					else if(geneAccordingToClinVar.contains(snpEffGene))
					{
						finalGeneSymbol = snpEffGene;
					}

					if(finalGeneSymbol != null)
					{
						if(index > 0)
						{
							variantToNonZeroSnpEffGeneIndex.put(chrPosRefAlt, index);
						}
						break;
					}
					index++;
				}
			}
			

			// still no joy?
			// there is a difference we cannot resolve nicely
			// so we concatenate the symbols from both
			if (finalGeneSymbol == null)
			{
				// simple case: we have 1 + 1
				if(genesAccordingToSnpEff.size() == 1)
				{
					finalGeneSymbol = geneAccordingToClinVar + "_" + genesAccordingToSnpEff.iterator().next();
				}
				
				
				//we make 1 exception for a mistake in clinvar: 'p.p.Arg801His' is actually gene 'RTEL1', but the notation is swapped
				//hopefully this will be fixed in the future so this exception is no longer needed (and removed)
				else if(geneAccordingToClinVar.equals("p.p.Arg801His"))
				{
					finalGeneSymbol = "RTEL1";
				}
				
				// if still not solved, we have 1 clinvar and 2+ snpeff symbols, so simply add them all to 1 big symbol
				// note 1: does not occur with clinvar november 2015 release
				// note 2: this would be slightly problematic because we don't know which SnpEff ANN field to consider later on!
				else
				{
					String snpEffGenes = "";
					for(String snpEffGene : genesAccordingToSnpEff)
					{
						snpEffGenes +=  "_" + snpEffGene;
					}
					finalGeneSymbol = geneAccordingToClinVar + snpEffGenes;
				}
			}		

			if (clinvarPatho.containsKey(finalGeneSymbol))
			{
				clinvarPatho.get(finalGeneSymbol).add(record);
			}
			else
			{
				List<Entity> variants = new ArrayList<Entity>();
				variants.add(record);
				clinvarPatho.put(finalGeneSymbol, variants);
			}
		}
		
		System.out.println("there are " + variantToNonZeroSnpEffGeneIndex.size() + " ClinVar variants for which the first SnpEff gene symbol is not the one matched to the ClinVar gene symbol");
		
		vcfRepo.close();
	}

	private void createMatchingExACsets(String exacLoc) throws Exception
	{
		Step4_Helper st4h = new Step4_Helper(variantToNonZeroSnpEffGeneIndex);
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
//			if(!gene.equals("TYROBP"))
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
			
			// include (biggest) part of exon the variant(s) are in, typical exon is 147 nt
			// http://nar.oxfordjournals.org/content/early/2012/07/11/nar.gks652.full
			leftMostPos = leftMostPos - 100;
			rightMostPos = rightMostPos + 100;

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
			
			//there are a lot of genes with only 1 pathogenic variant.. a bit silly to consider them seriously for calibration work
			//drop and report as N1
			if(clinvarPatho.get(gene).size() < 2)
			{
				droppedGenesClinVarTooFew++;
				//report exac impact and MAF (if overlaps with clinvar) anyway
				String exacImpact = exacVariants.size() > 0 ? "\t" + st4h.calculateImpactRatiosFromUnprocessedVariants(exacVariants).toString() : StringUtils.repeat("\t" + NA, 4);
				String maf = exacVariants.size() > 0 ? st4h.getExACMAFforUnprocessedClinvarVariant(clinvarPatho.get(gene).get(0), exacVariants) : NA;
				geneInfo.put(gene, "N1" + "\t" + clinvarPatho.get(gene).get(0).getString("#CHROM") + "\t" + clinvarPatho.get(gene).get(0).getString("POS") + "\t" + clinvarPatho.get(gene).get(0).getString("POS") + "\t" + exacVariants.size() + "\t" + clinvarPatho.get(gene).size() + "\t" + (maf.equals(NA) ? 0 : 1) + "\t" + 0 + "\t" + maf + exacImpact + "\t" + st4h.calculateImpactRatiosFromUnprocessedVariants(clinvarPatho.get(gene)) + StringUtils.repeat("\t" + NA, 4));
				continue;
			}

			System.out.println("\n#####\n");
			System.out.println(gene + " (" + index + " of " + clinvarPatho.keySet().size() + ") " + leftMostPos + " " + rightMostPos + " " + " has " + exacVariants.size());

			if (exacVariants.size() > 0)
			{
				//found out: which variants are only in ExAC, which only in ClinVar, which in both
				VariantIntersectResult vir = st4h.intersectVariants(exacVariants, clinvarPatho.get(gene));
				
				System.out.println("VariantIntersectResult for '"+gene+"', clinvaronly: " + vir.inClinVarOnly.size() + ", exaconly: " + vir.inExACOnly.size() + ", both: " + vir.inBoth_exac.size());
				
				//calculate MAF for shared variants, and use them to filter the other ExAC variants
				//this way, we use the overlap to determine a fair cutoff for the 'assumed benign' variants
				//if we have nothing to go on, we will set this to 0 and only select singleton variants
				double pathogenicMAF = st4h.calculatePathogenicMAF(vir.inBoth_exac, vir.inClinVarOnly.size());
				List<EntityPlus> exacFilteredByMAF = st4h.filterExACvariantsByMAF(vir.inExACOnly, pathogenicMAF);

				System.out.println("exaconly filtered down to " + exacFilteredByMAF.size() + " variants using pathogenic MAF " + pathogenicMAF);

				//calculate impact ratios over all clinvar variants, and use them to 'shape' the remaining ExAC variants
				//they must become a set that looks just like the ClinVar variants, including same distribution of impact types
				ImpactRatios pathoImpactRatio = st4h.calculateImpactRatios(Stream.concat(vir.inClinVarOnly.stream(), vir.inBoth_clinvar.stream()).collect(Collectors.toList()));
				String unfilteredExacImpactRatio = vir.inExACOnly.size() > 0 ? "\t" + st4h.calculateImpactRatios(vir.inExACOnly).toString() : StringUtils.repeat("\t" + NA, 4);
				
				if(exacFilteredByMAF.size() == 0)
				{
					geneInfo.put(gene, "T1" + "\t" + chrom + "\t" + leftMostPos + "\t" + rightMostPos + "\t" + exacVariants.size() + "\t" + clinvarPatho.get(gene).size() + "\t" + vir.inBoth_clinvar.size() + "\t" + 0 + "\t" + pathogenicMAF + unfilteredExacImpactRatio + "\t" + pathoImpactRatio.toString() + StringUtils.repeat("\t" + NA, 4));
					continue;
				}
				
				List<EntityPlus> exacFilteredByMAFandImpact = st4h.shapeExACvariantsByImpactRatios(exacFilteredByMAF, pathoImpactRatio);
				System.out.println("exaconly filtered down to " + exacFilteredByMAFandImpact.size() + " variants");
				
				
				if (exacFilteredByMAFandImpact.size() > 0)
				{
					passedGenes++;
					matchedExACvariants.put(gene, exacFilteredByMAFandImpact);
					matchedVariants += exacFilteredByMAFandImpact.size();
					//impacts AFTER impact correction and MAF filter has been applied
					ImpactRatios MAFandImpactFilteredExacImpactRatio = st4h.calculateImpactRatios(exacFilteredByMAFandImpact);
					geneInfo.put(gene, "Cx" + "\t" + chrom + "\t" + leftMostPos + "\t" + rightMostPos + "\t" + exacVariants.size() + "\t" + clinvarPatho.get(gene).size() + "\t" + vir.inBoth_clinvar.size() + "\t" + exacFilteredByMAFandImpact.size() + "\t" + pathogenicMAF + unfilteredExacImpactRatio + "\t" + pathoImpactRatio.toString().toString() + "\t" + MAFandImpactFilteredExacImpactRatio.toString());
				}
				else
				{
					//impacts BEFORE impact correction (which whould have set it to 0) but AFTER the MAF filter has been applied
					ImpactRatios MAFfilteredExacImpactRatio = st4h.calculateImpactRatios(exacFilteredByMAF);
					String cat = st4h.determineImpactFilterCat(MAFfilteredExacImpactRatio, pathoImpactRatio, pathogenicMAF);
					geneInfo.put(gene, cat + "\t" + chrom + "\t" + leftMostPos + "\t" + rightMostPos + "\t" + exacVariants.size() + "\t" + clinvarPatho.get(gene).size() + "\t" + vir.inBoth_clinvar.size() + "\t" + exacFilteredByMAF.size() + "\t" + pathogenicMAF + unfilteredExacImpactRatio + "\t" + pathoImpactRatio.toString() + "\t" + MAFfilteredExacImpactRatio.toString());
					droppedGenesNoMatchedVariants++;
				}
				
				

			}
			else
			{
				droppedGenesExACtooFew++;
				geneInfo.put(gene, "N2" + "\t" + chrom + "\t" + leftMostPos + "\t" + rightMostPos + "\t" + 0 + "\t" + clinvarPatho.get(gene).size() + "\t" + 0 + "\t" + 0 + "\t" + 0 + StringUtils.repeat("\t" + NA, 4)  + "\t" + st4h.calculateImpactRatiosFromUnprocessedVariants(clinvarPatho.get(gene)) + StringUtils.repeat("\t" + NA, 4));
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
		pw_geneInfo.println( "Gene" + "\t" + "Category" + "\t" + "Chr" + "\t" + "Start" + "\t" + "End" + "\t" + "NrOfPopulationVariants" + "\t" + "NrOfPathogenicVariants" + "\t" + "NrOfOverlappingVariants" + "\t" + "NrOfFilteredPopVariants" + "\t" + "PathoMAFThreshold" + "\t" + "PopImpactHighPerc" + "\t" + "PopImpactModeratePerc" + "\t" + "PopImpactLowPerc" + "\t" + "PopImpactModifierPerc" + "\t" + "PathoImpactHighPerc" + "\t" + "PathoImpactModeratePerc" + "\t" + "PathoImpactLowPerc" + "\t" + "PathoImpactModifierPerc" + "\t" + "PopImpactHighEq" + "\t" + "PopImpactModerateEq" + "\t" + "PopImpactLowEq" + "\t" + "PopImpactModifierEq");
		
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
