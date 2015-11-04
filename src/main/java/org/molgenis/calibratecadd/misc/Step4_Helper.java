package org.molgenis.calibratecadd.misc;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import org.apache.commons.math3.stat.descriptive.rank.Percentile;
import org.apache.commons.math3.stat.descriptive.rank.Percentile.EstimationType;
import org.molgenis.calibratecadd.structs.EntityPlus;
import org.molgenis.calibratecadd.structs.ImpactRatios;
import org.molgenis.calibratecadd.structs.VariantIntersectResult;
import org.molgenis.data.Entity;

public class Step4_Helper
{

	public static VariantIntersectResult intersectVariants(List<Entity> exacMultiAllelic, List<Entity> clinvar) throws Exception
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
			for(int altIndex = 0; altIndex < altSplit.length; altIndex++)
			{
				String alt = altSplit[altIndex];
				EntityPlus exacVariantCopy = new EntityPlus(exacVariant);
				exacVariantCopy.getKeyVal().put("ALT", alt);
				exacVariantCopy.getKeyVal().put("AF", Double.parseDouble(exacVariant.getString("AF").split(",", -1)[altIndex]));
				exac.add(exacVariantCopy);
			}
			//set original to null so we don't accidentally use it somewhere
			exacVariant.set("ALT", null);
			exacVariant.set("AF", null);
		}
		
//		System.out.println("expanded exac from " + exacMultiAllelic.size() + " to " + exac.size());
		
//		System.out.println("### BEFORE: ");
//		for(EntityPlus e : exac)
//		{
//			System.out.println(e.getE().getString("#CHROM") + "_" + e.getE().getString("POS") + "_" + e.getE().getString("REF") + "_" + e.getKeyVal().get("ALT").toString());
//		}
//		for(Entity e : clinvar)
//		{
//			System.out.println(e.getString("#CHROM") + "_" + e.getString("POS") + "_" + e.getString("REF") + "_" + e.get("ALT").toString());
//		}
		
		for (EntityPlus exacVariant : exac)
		{
			boolean exacVariantInClinVar = false;
			for (Entity clinvarVariant : clinvar)
			{
				// TODO
				// for now, we accept that we will miss *some* variants due to 2 reasons:
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
		
//		System.out.println("### AFTER: ");
//		for(EntityPlus e : inExAConly)
//		{
//			System.out.println(e.getE().getString("#CHROM") + "_" + e.getE().getString("POS") + "_" + e.getE().getString("REF") + "_" + e.getKeyVal().get("ALT"));
//		}
//		for(EntityPlus e : inClinVarOnly)
//		{
//			System.out.println(e.getE().getString("#CHROM") + "_" + e.getE().getString("POS") + "_" + e.getE().getString("REF") + "_" + e.getE().getString("ALT"));
//		}
//		for(EntityPlus e : inBoth_exac)
//		{
//			System.out.println(e.getE().getString("#CHROM") + "_" + e.getE().getString("POS") + "_" + e.getE().getString("REF") + "_" + e.getKeyVal().get("ALT"));
//		}
//		for(EntityPlus e : inBoth_clinvar)
//		{
//			System.out.println(e.getE().getString("#CHROM") + "_" + e.getE().getString("POS") + "_" + e.getE().getString("REF") + "_" + e.getE().getString("ALT"));
//		}
		
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
	
	
	
	public static double calculatePathogenicMAF(List<EntityPlus> exacVariants, int nrOfClinVarOnly)
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
		System.out.println("mafs: " + Arrays.toString(mafs));
		//R7 is the one used by R and Excel as default
		Percentile perc = new Percentile().withEstimationType(EstimationType.R_7);
		return perc.evaluate(mafs, 95);
	}



	public static List<EntityPlus> filterExACvariantsByMAF(List<EntityPlus> inExACOnly, double MAFthreshold) throws Exception
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
					
					//update the 'variant annotation' line 'CSQ' to match this alt
					//includes adding 'impact' for later use
					boolean success = updateCSQ(exacVariantPlus, alt, maf, AC_Adj);
					
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

	private static boolean updateCSQ(EntityPlus exacVariant, String altAllele, double maf, int AC_Adj) throws Exception
	{
		String csq = exacVariant.getE().getString("CSQ");	
//		System.out.println("LOOKING FOR ALT " + altAllele);
		
		boolean found = false;
		
		//multiple transcripts, with each multiple alleles
		for(String csqs : csq.split(",", -1))
		{
			String[] csqSplit = csqs.split("\\|", -1);
			
//			System.out.println("csqSplit[0]="+csqSplit[0]);
//			System.out.println("csqSplit[18]="+csqSplit[18]);
	
			if(csqSplit[0].equals(altAllele) && csqSplit[18].equals("YES"))
			{
				exacVariant.getE().set("CSQ", csqs);
				exacVariant.getE().set("ALT", altAllele);
				exacVariant.getE().set("AF", maf);
				exacVariant.getE().set("AC_Adj", AC_Adj);
				String impact = getHighestImpact(csqSplit[4]);
				exacVariant.getKeyVal().put(VEPimpactCategories.IMPACT, impact);
				found = true;
				break;
			}
		}
		
		if(!found)
		{
	//		System.out.println("could not return CSQ, no alt allele match for '"+altAllele+"' && 'YES' consensus for " + csq);
			return false;
		}
		else
		{
			return true;
		}
		
	}

	//e.g. "splice_acceptor_variant&non_coding_transcript_variant"
	private static String getHighestImpact(String csqConsequences) throws Exception
	{
		int highestImpactRank = -1;
		for(String consequence : csqConsequences.split("&", -1))
		{
			String impact = VEPimpactCategories.getImpact(consequence);
			int impactRank = VEPimpactCategories.getImpactRank(impact);
			if(impactRank > highestImpactRank)
			{
				highestImpactRank = impactRank;
			}
		}
		if(highestImpactRank == -1)
		{
			throw new Exception("no impact match on " + csqConsequences);
		}
		return highestImpactRank == 3 ? "HIGH" : highestImpactRank == 2 ? "MODERATE" : highestImpactRank == 1 ? "LOW" : "MODIFIER";
	}

	public static ImpactRatios calculateImpactRatios(List<EntityPlus> inClinVarOnly, List<EntityPlus> inBoth) throws Exception
	{
		int nrOfHigh = 0;
		int nrOfModerate = 0;
		int nrOfLow = 0;
		int nrOfModifier = 0;
		
		for(EntityPlus e : Stream.concat(inClinVarOnly.stream(), inBoth.stream()).collect(Collectors.toList()))
		{
			String impact = e.getE().getString("ANN").split("\\|")[2];
			if(impact.equals("HIGH"))
			{
				nrOfHigh++;
			}
			else if(impact.equals("MODERATE"))
			{
				nrOfModerate++;
			}
			else if(impact.equals("LOW"))
			{
				nrOfLow++;
			}
			else if(impact.equals("MODIFIER"))
			{
				nrOfModifier++;
			}
			else
			{
				throw new Exception("unrecognized impact: " + impact);
			}
		}
		
		double total = nrOfHigh + nrOfModerate + nrOfLow + nrOfModifier;
		double highPerc = nrOfHigh == 0 ? 0 :((double)nrOfHigh/total)*100.0;
		double modrPerc = nrOfModerate == 0 ? 0 : ((double)nrOfModerate/total)*100.0;
		double lowPerc = nrOfLow == 0 ? 0 : ((double)nrOfLow/total)*100.0;
		double modfPerc = nrOfModifier == 0 ? 0 : ((double)nrOfModifier/total)*100;
		
		ImpactRatios ir = new ImpactRatios(highPerc, modrPerc, lowPerc, modfPerc);
		
	//	System.out.println("counts: high=" + nrOfHigh + ", modr=" + nrOfModerate + ", low=" + nrOfLow + ", modf=" + nrOfModifier);

		return ir;
	}

	/**
	 * TODO: this function has an interesting side effect: when there are (only) HIGH effect variants in clinvar, but only MODERATE (or LOW/MODF) variants in ExAC, the matching fails..
	 * However, we do learn that apparently a HIGH impact variant is pathogenic, whereas non-HIGH are tolerated to some point. Even though we cannot calibrate CADD, this knowledge is
	 * just as useful and we should capture and report it :)
	 */
	public static List<EntityPlus> shapeExACvariantsByImpactRatios(List<EntityPlus> exacFilteredByMAF, ImpactRatios ir) throws Exception
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
			String impact = e.getKeyVal().get(VEPimpactCategories.IMPACT).toString();
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
	
	
	private static List<EntityPlus> scaledownVariantList(List<EntityPlus> variants, int amountToRemove) throws Exception
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

	public static String determineImpactFilterCat(List<EntityPlus> exacFilteredByMAF, ImpactRatios ir, double pathoMAF) throws Exception
	{
		
		//first, just count the impact categories like we do for clinvar
		int nrOfHighInExAC = 0;
		int nrOfModerateInExAC = 0;
		int nrOfLowInExAC = 0;
		int nrOfModifierInExAC = 0;
		
		for(EntityPlus e : exacFilteredByMAF)
		{
			String impact = e.getKeyVal().get(VEPimpactCategories.IMPACT).toString();
			if(impact.equals("HIGH"))
			{
				nrOfHighInExAC++;
			}
			else if(impact.equals("MODERATE"))
			{
				nrOfModerateInExAC++;
			}
			else if(impact.equals("LOW"))
			{
				nrOfLowInExAC++;
			}
			else if(impact.equals("MODIFIER"))
			{
				nrOfModifierInExAC++;
			}
			else
			{
				throw new Exception("unrecognized impact: " + impact);
			}
		}

		String details = "\t" + "pathogenic MAF: " + pathoMAF + ", pathogenic impact ratio: " + ir.toString() + ", exac impact counts: high="+nrOfHighInExAC+", modr="+nrOfModerateInExAC+", low="+nrOfLowInExAC + ", modf="+nrOfModifierInExAC;
		
		if(ir.high > 0 && nrOfHighInExAC == 0)
		{
			return "I1" + details;
		}
		else if(ir.moderate > 0 && nrOfHighInExAC == 0 && nrOfModerateInExAC == 0)
		{
			return "I2" + details;
		}
		else if(ir.low > 0 && nrOfHighInExAC == 0 && nrOfModerateInExAC == 0 && nrOfLowInExAC == 0)
		{
			return "I3" + details;
		}
		else
		{
			return "T2" + details;
		}

	}
	
}
