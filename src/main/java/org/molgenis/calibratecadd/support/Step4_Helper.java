package org.molgenis.calibratecadd.support;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import org.apache.commons.math3.stat.descriptive.rank.Percentile;
import org.apache.commons.math3.stat.descriptive.rank.Percentile.EstimationType;
import org.molgenis.data.Entity;

public class Step4_Helper
{
	private HashMap<String, Integer> variantToNonZeroSnpEffGeneIndex;
	
	/**
	 * Constructor
	 * @param variantToNonZeroSnpEffGeneIndex
	 */
	public Step4_Helper(HashMap<String, Integer> variantToNonZeroSnpEffGeneIndex)
	{
		this.variantToNonZeroSnpEffGeneIndex = variantToNonZeroSnpEffGeneIndex;
	}

	/**
	 * If these ExAC and ClinVar variants have a matching ALT allele, report the ExAC AF, else return 0
	 * @param clinvarVariant
	 * @return
	 * @throws Exception
	 */
	public String getExACMAFforUnprocessedClinvarVariant(Entity clinvarVariant, List<Entity> exacVariants) throws Exception
	{
		for(Entity exacVar : exacVariants)
		{
			if (exacVar.getString("#CHROM").equals(clinvarVariant.getString("#CHROM"))
					&& exacVar.getString("POS").equals(clinvarVariant.getString("POS"))
					&& exacVar.getString("REF").equals(clinvarVariant.getString("REF"))
					)
			{
				String clinvarAlt = clinvarVariant.getString("ALT");
				String[] altSplit = exacVar.getString("ALT").split(",", -1);
				for(int altIndex = 0; altIndex < altSplit.length; altIndex++)
				{
					if(clinvarAlt.equals(altSplit[altIndex]))
					{
						return exacVar.getString("AF").split(",", -1)[altIndex];
					}
				}
			}
		}
		return "0";
	}
	
	/**
	 * Take 2 lists and return 4: exac-only, clinvar-only, exac-inboth, clinvar-inboth.
	 * @param exacMultiAllelic
	 * @param clinvar
	 * @return
	 * @throws Exception
	 */
	public VariantIntersectResult intersectVariants(List<Entity> exacMultiAllelic, List<Entity> clinvar) throws Exception
	{
		List<EntityPlus> inExAConly = new ArrayList<EntityPlus>();
		List<EntityPlus> inClinVarOnly = new ArrayList<EntityPlus>();
		List<EntityPlus> inBoth_exac = new ArrayList<EntityPlus>();
		List<EntityPlus> inBoth_clinvar = new ArrayList<EntityPlus>();
		
		//preprocess: expand multiallelic ExaC variants into seperate variants, for 'easy of looping'
		//update the alt allele and AF. This only goes for ExAC variants because they can be multi-allelic, whereas ClinVar variants are not.
		//example of where this goes wrong if we don't do this: 6:32007887 . Here, there is an ExAC variant G -> T,A and ClinVar G -> C and G -> T.
		//move ALT and AF fields to keyVal map
		List<EntityPlus> exac = new ArrayList<EntityPlus>();
		for (Entity exacVariant : exacMultiAllelic)
		{
			String[] altSplit = exacVariant.getString("ALT").split(",", -1);
			Set<String> altsSeenForVariant = new HashSet<String>();
			for(int altIndex = 0; altIndex < altSplit.length; altIndex++)
			{
				String alt = altSplit[altIndex];
				EntityPlus exacVariantCopy = new EntityPlus(exacVariant);
				//sanity check
				if(altsSeenForVariant.contains(alt))
				{
					throw new Exception("Same alt seen twice for " + exacVariant.toString());
				}
				altsSeenForVariant.add(alt);
				exacVariantCopy.getKeyVal().put("ALT", alt);
				exacVariantCopy.getKeyVal().put("AF", Double.parseDouble(exacVariant.getString("AF").split(",", -1)[altIndex]));
				exac.add(exacVariantCopy);
			}
			//set original to null so we don't accidentally use it somewhere
			exacVariant.set("ALT", null);
			exacVariant.set("AF", null);
		}
		
		for (EntityPlus exacVariant : exac)
		{
			boolean exacVariantInClinVar = false;
			for (Entity clinvarVariant : clinvar)
			{
				// TODO
				// for now, we accept that we might miss *some* variants due to 2 reasons:
				// 1) offset positions due to complex indels
				// 2) alternative notation of indels, e.g.: consider this variant: 1 6529182 . TTCCTCC TTCC
				// you will find that it is seen in ExAC: 1 6529182 . TTCCTCCTCC TTCCTCC,TTCC,T,TTCCTCCTCCTCC,TTCCTCCTCCTCCTCC,TTCCTCCTCCTCCTCCTCCTCC
				// but there denoted as "TTCCTCCTCC/TTCCTCC"...
				if (exacVariant.getE().getString("#CHROM").equals(clinvarVariant.getString("#CHROM"))
						&& exacVariant.getE().getString("POS").equals(clinvarVariant.getString("POS"))
						&& exacVariant.getE().getString("REF").equals(clinvarVariant.getString("REF"))
						&& exacVariant.getKeyVal().get("ALT").toString().equals(clinvarVariant.getString("ALT")))
				{
					inBoth_exac.add(exacVariant);
					inBoth_clinvar.add(new EntityPlus(clinvarVariant));
					exacVariantInClinVar = true;
				}
			}
			if(!exacVariantInClinVar)
			{
				inExAConly.add(exacVariant);
			}
		}
		
		// now have have the list of variants that are shared
		// do a pass of clinvar variants and find out which are not shared
		for (Entity clinVarVariant : clinvar)
		{
			if(!EntityPlus.contains(inBoth_clinvar, clinVarVariant))
			{
				inClinVarOnly.add(new EntityPlus(clinVarVariant));
			}
		}
				
		//sanity checks
		if(inBoth_clinvar.size() != inBoth_exac.size())
		{
			throw new Exception("inBoth sizes not equal: " + inBoth_clinvar.size() + " vs " + inBoth_exac.size());
		}
		if(exac.size()+clinvar.size() != inExAConly.size()+inClinVarOnly.size()+inBoth_exac.size()+inBoth_clinvar.size())
		{
			throw new Exception("Sizes dont add up: " + exac.size() + "+" + clinvar.size() + " != " + inExAConly.size() + "+" +inClinVarOnly.size() + "+" + inBoth_exac.size() + "+" + inBoth_clinvar.size());
		}
		
		return new VariantIntersectResult(inBoth_exac, inBoth_clinvar, inClinVarOnly, inExAConly);
	}
	
	/**
	 * Given a set of ExAC variants and a number of variants that was unique to ClinVar, calculate 95th percentile MAF
	 * @param exacVariants
	 * @param nrOfClinVarOnly
	 * @return
	 */
	public double calculatePathogenicMAF(List<EntityPlus> exacVariants, int nrOfClinVarOnly)
	{
		if(exacVariants.size() == 0)
		{
			return 0.0;
		}
		double[] mafs = new double[exacVariants.size() + nrOfClinVarOnly];
		for(int i = 0; i < exacVariants.size(); i++)
		{
			mafs[i] = (double) exacVariants.get(i).getKeyVal().get("AF");
		}
		for(int i = exacVariants.size(); i < exacVariants.size() + nrOfClinVarOnly; i++)
		{
			mafs[i] = 0.0;
		}
		//R7 is the one used by R and Excel as default
		Percentile perc = new Percentile().withEstimationType(EstimationType.R_7);
		return perc.evaluate(mafs, 95);
	}


	/**
	 * Given a MAF threshold, filter a set of ExAC variants
	 * If MAF = 0, we only select singleton variants (AC_Adj = 1)
	 * @param inExACOnly
	 * @param MAFthreshold
	 * @return
	 * @throws Exception
	 */
	public List<EntityPlus> filterExACvariantsByMAF(List<EntityPlus> inExACOnly, double MAFthreshold, String gene) throws Exception
	{
		List<EntityPlus> res = new ArrayList<EntityPlus>();
		
		filterVariants:
		for(EntityPlus exacVariantPlus : inExACOnly)
		{
	
			String[] altSplit = exacVariantPlus.getKeyVal().get("ALT").toString().split(",", -1);
			for(int altIndex = 0; altIndex < altSplit.length; altIndex++)
			{
				String alt = altSplit[altIndex];
				//System.out.println("AF field: " + exacVariantPlus.getE().getString("AF"));
			
				double maf = Double.parseDouble(exacVariantPlus.getKeyVal().get("AF").toString().split(",",-1)[altIndex]);				
				int AC_Adj = Integer.parseInt(exacVariantPlus.getE().getString("AC_Adj").split(",",-1)[altIndex]);
		
				//we consider each alt allele as a possible 'keep' or 'ditch'
				//though we only keep 1 alt allele if that one is a match
				boolean keep = false;
				
				//the clinvar variants were all 'singletons', so we only select singletons from exac
				if(MAFthreshold == 0.0 && AC_Adj == 1)
				{
					keep = true;
				}
				//else it must be under/equal to MAF threshold
				else if(maf <= MAFthreshold)
				{
					keep = true;
				}
				
				//if keep: we have to update this variant to remove any 'ditched' alternative alleles!
				if(keep)
				{
					
					//update the 'variant annotation' line 'ANN' to match this alt
					//includes adding 'impact' for later use
					boolean success = updateANN(exacVariantPlus, alt, maf, AC_Adj, gene);
					
					if(success)
					{
						res.add(exacVariantPlus);

						//don't consider other alt alleles, just stick with the one we found first and continue
						continue filterVariants;
					}
				}
			}
		}
//		System.out.println("RETURNING " + res.get(0).getKeyVal().get(VEPimpactCategories.IMPACT).toString());
		return res;
	}

	/**
	 * Helper function to update ANN (annotation) field to only contain the ANN for the matching
	 * alt allele, containing the highest impact effect, and for canonical transcripts only
	 * @param exacVariant
	 * @param altAllele
	 * @param maf
	 * @param AC_Adj
	 * @return
	 * @throws Exception
	 */
	private boolean updateANN(EntityPlus exacVariant, String altAllele, double maf, int AC_Adj, String gene) throws Exception
	{
		String alt = exacVariant.getKeyVal().get("ALT").toString();
		String fullAnn = exacVariant.getE().getString("ANN");
//		System.out.println("FULLANN " + fullAnn);
		String ann = GavinUtils.getAnn(fullAnn, gene, alt);
		if(ann == null)
		{
			return false;
		}
		exacVariant.getE().set("ANN", fullAnn); //dont put alt-specific ANN in, because we iterate over multiple alts..
		exacVariant.getE().set("ALT", altAllele);
		exacVariant.getE().set("AF", maf);
		exacVariant.getE().set("AC_Adj", AC_Adj);
	//	System.out.println("ANN = " + ann);
		String impact = GavinUtils.getImpact(fullAnn, gene, alt).toString();
		exacVariant.getKeyVal().put("IMPACT", impact);

		return true;

	}

	public ImpactRatios calculateImpactRatiosFromUnprocessedVariants(List<Entity> entities) throws Exception
	{
		//liftover entity to entityplus and count
		List<EntityPlus> clinvarVariants = new ArrayList<EntityPlus>();
		for(Entity e : entities)
		{
			clinvarVariants.add(new EntityPlus(e));
		}
		return calculateImpactRatios(clinvarVariants);
	}
	
	/**
	 * Pass existing counts by reference, and increment one value based on the impact:
	 * counts[0] = "HIGH"
	 * counts[1] = "MODERATE"
	 * counts[2] = "LOW"
	 * counts[3] = "MODIFIER"
	 * @param counts
	 * @param impact
	 * @throws Exception
	 *
	 */
	public void countImpacts(Integer[] counts, String impact) throws Exception
	{
		if(impact.equals("HIGH"))
		{
			counts[0]++;
		}
		else if(impact.equals("MODERATE"))
		{
			counts[1]++;
		}
		else if(impact.equals("LOW"))
		{
			counts[2]++;
		}
		else if(impact.equals("MODIFIER"))
		{
			counts[3]++;
		}
		else
		{
			throw new Exception("unrecognized impact: " + impact);
		}
	}
	
	/**
	 * Determine ratio of impacts for a list of variants
	 * @param entities
	 * @return
	 * @throws Exception
	 */
	public ImpactRatios calculateImpactRatios(List<EntityPlus> entities) throws Exception
	{
		Integer[] impactCounts = new Integer[]{ 0, 0, 0, 0};
		
		for(EntityPlus e : entities)
		{
			
			if(e.getE().getString("ANN") != null)
			{
				String impact = null;
				//clinvar
				if(e.getE().getString("ANN") != null)
				{
					String[] multiAnn = e.getE().getString("ANN").split(",", -1);
					
					if(multiAnn.length > 1)
					{
						//special: the ANN field for the gene symbol we mapped is not present on index 0
						//get it from the list of 'exceptions' we stored earlier
						String chrPosRefAlt = e.getE().getString("#CHROM") + "_" + e.getE().getString("POS") + "_" + e.getE().getString("REF") + "_" + e.getE().getString("ALT");
						if(variantToNonZeroSnpEffGeneIndex.containsKey(chrPosRefAlt))
						{
							impact = multiAnn[variantToNonZeroSnpEffGeneIndex.get(chrPosRefAlt)].split("\\|", -1)[2];
						}
						else
						{
							impact = multiAnn[0].split("\\|", -1)[2];
						}
					}
					else
					{
						impact =  multiAnn[0].split("\\|", -1)[2];
					}
					
				}
				else
				{
					throw new Exception("should not be reached");
				}
				
				countImpacts(impactCounts, impact);
				
			}
			else
			{
				throw new Exception("should no longer be reached when ExAC is annotated by SnpEff");
			}
		}
		
		double total = impactCounts[0] + impactCounts[1] + impactCounts[2] + impactCounts[3];
		double highPerc = impactCounts[0] == 0 ? 0 :((double)impactCounts[0]/total)*100.0;
		double modrPerc = impactCounts[1] == 0 ? 0 : ((double)impactCounts[1]/total)*100.0;
		double lowPerc = impactCounts[2] == 0 ? 0 : ((double)impactCounts[2]/total)*100.0;
		double modfPerc = impactCounts[3] == 0 ? 0 : ((double)impactCounts[3]/total)*100;
		
		ImpactRatios ir = new ImpactRatios(highPerc, modrPerc, lowPerc, modfPerc);
		
	//	System.out.println("counts: high=" + nrOfHigh + ", modr=" + nrOfModerate + ", low=" + nrOfLow + ", modf=" + nrOfModifier);

		return ir;
	}

	/**
	 * TODO: this function has an interesting side effect: when there are (only) HIGH effect variants in clinvar, but only MODERATE (or LOW/MODF) variants in ExAC, the matching fails..
	 * However, we do learn that apparently a HIGH impact variant is pathogenic, whereas non-HIGH are tolerated to some point. Even though we cannot calibrate CADD, this knowledge is
	 * just as useful and we should capture and report it :)
	 */
	public List<EntityPlus> shapeExACvariantsByImpactRatios(List<EntityPlus> exacFilteredByMAF, ImpactRatios ir) throws Exception
	{
		List<EntityPlus> highImpactVariants = new ArrayList<EntityPlus>();
		List<EntityPlus> modrImpactVariants = new ArrayList<EntityPlus>();
		List<EntityPlus> lowImpactVariants = new ArrayList<EntityPlus>();
		List<EntityPlus> modfImpactVariants = new ArrayList<EntityPlus>();
		
		//first, just count the impact categories like we do for clinvar
		int nrOfHigh = 0;
		int nrOfModerate = 0;
		int nrOfLow = 0;
		int nrOfModifier = 0;
		for(EntityPlus e : exacFilteredByMAF)
		{
	//		System.out.println(e.getKeyVal().toString());
			String impact = e.getKeyVal().get("IMPACT").toString();
			if(impact.equals("HIGH"))
			{
				highImpactVariants.add(e);
				nrOfHigh++;
			}
			else if(impact.equals("MODERATE"))
			{
				modrImpactVariants.add(e);
				nrOfModerate++;
			}
			else if(impact.equals("LOW"))
			{
				lowImpactVariants.add(e);
				nrOfLow++;
			}
			else if(impact.equals("MODIFIER"))
			{
				modfImpactVariants.add(e);
				nrOfModifier++;
			}
			else
			{
				throw new Exception("unrecognized impact: " + impact);
			}
		}
		
		System.out.println("counting exac impacts: high="+nrOfHigh+", modr="+nrOfModerate+", low="+nrOfLow + ", modf="+nrOfModifier);
		
		//tackle:
		//we have impact ratios, e.g.: [high=40, moderate=53, low=7, modifier=0]
		//counting exac impacts: high=25, modr=230, low=144, modf=235
		//limiting impact type: high with 25, total set we want: 25*(100/40) = 62.5
		//fill the rest: 62.5*(53/100) = 33.125 moderate impact ones, 4.375 low impact, 0 modifier
		//33+4+25 = 62, which is fine, just round up/down each impact type
		
		//so! little bit tricky: we must get the difference between 'initial' vs 'scaled' for 3 categories, using 1 as 'scaling reference'
		//if all 3 are negative, it means this is the correct 'scaling reference' because we can remove variants but not add variants!
		
		//alright.. let's test if 'high' is our scaling reference:
		int highScaleModrDiff = -1, highScaleLowDiff = -1, highScaleModfDiff = -1;
		if(ir.high != 0)
		{
			highScaleModrDiff = (int)Math.round(nrOfModerate-(nrOfHigh*(ir.moderate/ir.high)));
			highScaleLowDiff = (int)Math.round(nrOfLow-(nrOfHigh*(ir.low/ir.high)));
			highScaleModfDiff = (int)Math.round(nrOfModifier-(nrOfHigh*(ir.modifier/ir.high)));
		}
		
		
		//no? then check if we should scale on 'moderate'
		int modrScaleHighDiff = -1, modrScaleLowDiff = -1, modrScaleModfDiff = -1;
		if(ir.moderate != 0)
		{
			modrScaleHighDiff = (int)Math.round(nrOfHigh-(nrOfModerate*(ir.high/ir.moderate)));
			modrScaleLowDiff = (int)Math.round(nrOfLow-(nrOfModerate*(ir.low/ir.moderate)));
			modrScaleModfDiff = (int)Math.round(nrOfModifier-(nrOfModerate*(ir.modifier/ir.moderate)));
		}
		
		//no? then check if we should scale on 'low'
		int lowScaleHighDiff = -1, lowScaleModrDiff = -1, lowScaleModfDiff = -1;
		if(ir.low != 0)
		{
			lowScaleHighDiff = (int)Math.round(nrOfHigh-(nrOfLow*(ir.high/ir.low)));
			lowScaleModrDiff = (int)Math.round(nrOfModerate-(nrOfLow*(ir.moderate/ir.low)));
			lowScaleModfDiff = (int)Math.round(nrOfModifier-(nrOfLow*(ir.modifier/ir.low)));
		}
		
		//no? then check if we should scale on 'modifier'
		int modfScaleHighDiff = -1, modfScaleModrDiff = -1, modfScaleLowDiff = -1;
		if(ir.modifier != 0)
		{
			modfScaleHighDiff = (int)Math.round(nrOfHigh-(nrOfModifier*(ir.high/ir.modifier)));
			modfScaleModrDiff = (int)Math.round(nrOfModerate-(nrOfModifier*(ir.moderate/ir.modifier)));
			modfScaleLowDiff = (int)Math.round(nrOfLow-(nrOfModifier*(ir.low/ir.modifier)));
		}
		
		System.out.println("scaling subtractions for HIGH: moderate=" + highScaleModrDiff + ", low=" + highScaleLowDiff + ", modifier=" + highScaleModfDiff);
		System.out.println("scaling subtractions for MODERATE: high=" + modrScaleHighDiff + ", low=" + modrScaleLowDiff + ", modifier=" + modrScaleModfDiff);
		System.out.println("scaling subtractions for LOW: high=" + lowScaleHighDiff + ", moderate=" + lowScaleModrDiff + ", modifier=" + lowScaleModfDiff);
		System.out.println("scaling subtractions for MODIFIER: high=" + modfScaleHighDiff + ", moderate=" + modfScaleModrDiff + ", low=" + modfScaleLowDiff);
		
		int removeFromHigh = 0, removeFromModerate = 0, removeFromLow = 0, removeFromModifier = 0;
		
		//multiple solutions? pick one with LEAST amount of deleted elements
		//this happens when a category has 0 variants, the other categories all get scaled to 0 as well.. (e.g. from 10% to 50% of 0 is still 0)
		//however, this will result in a much bigger loss because all variants get deleted this way
		int leastLossSoFar = -1;
		
		if(highScaleModrDiff >= 0 && highScaleLowDiff >= 0 && highScaleModfDiff >= 0)
		{
		//	System.out.println("we must scale on HIGH impact using " + highScaleModrDiff + ", " + highScaleLowDiff + ", " + highScaleModfDiff);	
			int loss = highScaleModrDiff + highScaleLowDiff + highScaleModfDiff;
			if(leastLossSoFar == -1 || loss < leastLossSoFar)
			{
				leastLossSoFar = loss;
				removeFromHigh = 0;
				removeFromModerate = highScaleModrDiff;
				removeFromLow = highScaleLowDiff;
				removeFromModifier = highScaleModfDiff;
				System.out.println("scaling on HIGH is an option, with loss = " + loss);
			}
		}
		
		if(modrScaleHighDiff >= 0 && modrScaleLowDiff >= 0 && modrScaleModfDiff >= 0)
		{
			int loss = modrScaleHighDiff + modrScaleLowDiff + modrScaleModfDiff;
			if(leastLossSoFar == -1 || loss < leastLossSoFar)
			{
				leastLossSoFar = loss;
				removeFromHigh = modrScaleHighDiff;
				removeFromModerate = 0;
				removeFromLow = modrScaleLowDiff;
				removeFromModifier = modrScaleModfDiff;
				System.out.println("scaling on MODERATE is a (better) option, with loss = " + loss);
			}
//			System.out.println("we must scale on MODERATE impact using " + modrScaleHighDiff + ", " + modrScaleLowDiff + ", " + modrScaleModfDiff);
			
		}
		if(lowScaleHighDiff >= 0 && lowScaleModrDiff >= 0 && lowScaleModfDiff >= 0)
		{
			int loss = lowScaleHighDiff + lowScaleModrDiff + lowScaleModfDiff;
			if(leastLossSoFar == -1 || loss < leastLossSoFar)
			{
				leastLossSoFar = loss;
				removeFromHigh = lowScaleHighDiff;
				removeFromModerate = lowScaleModrDiff;
				removeFromLow = 0;
				removeFromModifier = lowScaleModfDiff;
				System.out.println("scaling on LOW is a (better) option, with loss = " + loss);
			}
			
		}
		if(modfScaleHighDiff >= 0 && modfScaleModrDiff >= 0 && modfScaleLowDiff >= 0)
		{
			int loss = modfScaleHighDiff + modfScaleModrDiff + modfScaleLowDiff;
			if(leastLossSoFar == -1 || loss < leastLossSoFar)
			{
				leastLossSoFar = loss;
				removeFromHigh = modfScaleHighDiff;
				removeFromModerate = modfScaleModrDiff;
				removeFromLow = modfScaleLowDiff;
				removeFromModifier = 0;
				System.out.println("scaling on MODIFIER is a (better) option, with loss = " + loss);
			}
		//	System.out.println("we must scale on MODIFIER impact using " + modfScaleHighDiff + ", " + modfScaleModrDiff + ", " + modfScaleLowDiff);
			
		}
		if(leastLossSoFar == -1)
		{
			throw new Exception("could not figure out scaling!");
		}
		
		
		System.out.println("removing from high: " + removeFromHigh + ", moderate: " + removeFromModerate + ", low: " + removeFromLow + ", modf: " + removeFromModifier);
		
		List<EntityPlus> highScaledDown = scaledownVariantList(highImpactVariants, removeFromHigh);
		List<EntityPlus> modrScaledDown = scaledownVariantList(modrImpactVariants, removeFromModerate);
		List<EntityPlus> lowScaledDown = scaledownVariantList(lowImpactVariants, removeFromLow);
		List<EntityPlus> modfScaledDown = scaledownVariantList(modfImpactVariants, removeFromModifier);
		
		List<EntityPlus> scaledDownExACvariants = new ArrayList<EntityPlus>();
		scaledDownExACvariants.addAll(highScaledDown);
		scaledDownExACvariants.addAll(modrScaledDown);
		scaledDownExACvariants.addAll(lowScaledDown);
		scaledDownExACvariants.addAll(modfScaledDown);
		
		return scaledDownExACvariants;
	}
	
	
	private List<EntityPlus> scaledownVariantList(List<EntityPlus> variants, int amountToRemove) throws Exception
	{
		//remove nothing
		if(amountToRemove == 0)
		{
			return variants;
		}
		
		//remove all
		int size = variants.size();
		if(size-amountToRemove == 0)
		{
			return new ArrayList<EntityPlus>();
		}

		//'how often does the final size fit within the list of variants? e.g. want 20 out of 190 variants = 9x
		//this means we will step through the variant list in steps of 9, to get 'even coverage'
		int div = Math.floorDiv(size, size-amountToRemove);
		List<EntityPlus> res = new ArrayList<EntityPlus>();
		for(int step = 0; step < size; step += div)
		{
			res.add(variants.get(step));
		}
		
	//	System.out.println("DIV: " + div);
		
		if(res.size() < size-amountToRemove)
		{
			throw new Exception("result too few! need" + (size-amountToRemove) + " variants, got " + res.size());
		}
		
		return res;
	}

	public String determineImpactFilterCat(ImpactRatios exacImpactRatio, ImpactRatios pathoImpactRatio, double pathoMAF) throws Exception
	{
		
		if(pathoImpactRatio.high > 0 && exacImpactRatio.high == 0)
		{
			return "I1";
		}
		else if(pathoImpactRatio.moderate > 0 && exacImpactRatio.high == 0 && exacImpactRatio.moderate == 0)
		{
			return "I2";
		}
		else if(pathoImpactRatio.low > 0 && exacImpactRatio.high == 0 && exacImpactRatio.moderate == 0 && exacImpactRatio.low == 0)
		{
			return "I3";
		}
		else
		{
			return "T2";
		}

	}
	
}
