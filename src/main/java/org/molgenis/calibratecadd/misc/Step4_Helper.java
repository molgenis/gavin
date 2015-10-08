package org.molgenis.calibratecadd.misc;

import java.util.ArrayList;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import org.apache.commons.math3.stat.descriptive.rank.Median;
import org.molgenis.calibratecadd.structs.EntityPlus;
import org.molgenis.calibratecadd.structs.ImpactRatios;
import org.molgenis.calibratecadd.structs.VariantIntersectResult;
import org.molgenis.data.Entity;
import org.molgenis.data.support.DefaultEntityMetaData;

public class Step4_Helper
{

	public static VariantIntersectResult intersectVariants(List<Entity> exac, List<Entity> clinvar) throws Exception
	{
		List<Entity> inExAConly = new ArrayList<Entity>();
		List<Entity> inClinVarOnly = new ArrayList<Entity>();
		List<Entity> inBoth_exac = new ArrayList<Entity>();
		List<Entity> inBoth_clinvar = new ArrayList<Entity>();
		
		for (Entity exacVariant : exac)
		{
			boolean exacVariantInClinVar = false;
			
			checkClinVarVariants:
			for (Entity clinvarVariant : clinvar)
			{

				// TODO
				// for now, we accept that we will miss *some* variants due to 2 reasons:
				// 1) offset positions due to complex indels
				// 2) alternative notation of indels, e.g.: consider this variant: 1 6529182 . TTCCTCC TTCC
				// you will find that it is seen in ExAC: 1 6529182 . TTCCTCCTCC TTCCTCC,TTCC,T,TTCCTCCTCCTCC,TTCCTCCTCCTCCTCC,TTCCTCCTCCTCCTCCTCCTCC
				// but there denoted as "TTCCTCCTCC/TTCCTCC"...
				if (exacVariant.getString("#CHROM").equals(clinvarVariant.getString("#CHROM"))
						&& exacVariant.getString("POS").equals(clinvarVariant.getString("POS")) && exacVariant.getString("REF").equals(clinvarVariant.getString("REF")))
				{
					String[] altSplit = exacVariant.getString("ALT").split(",", -1);
					for(int altIndex = 0; altIndex < altSplit.length; altIndex++)
					{
						String alt = altSplit[altIndex];
						if (alt.equals(clinvarVariant.getString("ALT")))
						{
//							System.out.println("MATCH ON chrom/pos/ref/alt:");
//							System.out.println(exacVariant.toString());
//							System.out.println(clinvarVariant.toString());
							
							//edit the variant if we have >1 alt so we only keep this particular alt allele and matching AF field
							if(altSplit.length > 1)
							{
								exacVariant.set("ALT", alt);
								exacVariant.set("AF", exacVariant.getString("AF").split(",", -1)[altIndex]);
//								System.out.println("CHANGED TO: ");
//								System.out.println(exacVariant.toString());
							}
							
							// TODO also merge data instead of 2 lists?
							inBoth_exac.add(exacVariant);
							inBoth_clinvar.add(clinvarVariant);
							exacVariantInClinVar = true;
							break checkClinVarVariants;
						}
						
					}
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
			if(!inBoth_clinvar.contains(clinVarVariant))
			{
				inClinVarOnly.add(clinVarVariant);
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
	
	
	
	public static double calculateMedianMAF(List<Entity> exacVariants)
	{
		if(exacVariants.size() == 0)
		{
			return 0.0;
		}
		double[] mafs = new double[exacVariants.size()];
		for(int i = 0; i < exacVariants.size(); i++)
		{
			mafs[i] = exacVariants.get(i).getDouble("AF");
		}
		Median m = new Median();
		return m.evaluate(mafs);
	}



	public static List<EntityPlus> filterExACvariantsByMAF(List<Entity> inExACOnly, double MAFthreshold) throws Exception
	{
		List<EntityPlus> res = new ArrayList<EntityPlus>();
		
		filterVariants:
		for(Entity exacVariant : inExACOnly)
		{
			EntityPlus exacVariantPlus = new EntityPlus(exacVariant);
			
			String[] altSplit = exacVariantPlus.getE().getString("ALT").split(",", -1);
			for(int altIndex = 0; altIndex < altSplit.length; altIndex++)
			{
				String alt = altSplit[altIndex];
				//System.out.println("AF field: " + exacVariantPlus.getE().getString("AF"));
			
				double maf = Double.parseDouble(exacVariantPlus.getE().getString("AF").split(",",-1)[altIndex]);				
				int AC_Adj = Integer.parseInt(exacVariantPlus.getE().getString("AC_Adj").split(",",-1)[altIndex]);
		
				//we consider each alt allele as a possible 'keep' or 'ditch'
				//though we only keep 1 alt allele if that one is a match
				boolean keep = false;
				
				//the clinvar variants were all 'singletons', so we only select singletons from exac
				if(maf == 0.0 && AC_Adj == 1)
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

	public static boolean updateCSQ(EntityPlus exacVariant, String altAllele, double maf, int AC_Adj) throws Exception
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
	public static String getHighestImpact(String csqConsequences) throws Exception
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

	public static ImpactRatios calculateImpactRatios(List<Entity> inClinVarOnly, List<Entity> inBoth) throws Exception
	{
		int nrOfHigh = 0;
		int nrOfModerate = 0;
		int nrOfLow = 0;
		int nrOfModifier = 0;
		
		for(Entity e : Stream.concat(inClinVarOnly.stream(), inBoth.stream()).collect(Collectors.toList()))
		{
			String impact = e.getString("ANN").split("\\|")[2];
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
		int highPerc = nrOfHigh == 0 ? 0 : (int)Math.round(((double)nrOfHigh/total)*100.0);
		int modrPerc = nrOfModerate == 0 ? 0 : (int)Math.round(((double)nrOfModerate/total)*100.0);
		int lowPerc = nrOfLow == 0 ? 0 : (int)Math.round(((double)nrOfLow/total)*100.0);
		int modfPerc = nrOfModifier == 0 ? 0 : (int)Math.round(((double)nrOfModifier/total)*100);
		
		ImpactRatios ir = new ImpactRatios(highPerc, modrPerc, lowPerc, modfPerc);
		
	//	System.out.println("counts: high=" + nrOfHigh + ", modr=" + nrOfModerate + ", low=" + nrOfLow + ", modf=" + nrOfModifier);

		return ir;
	}



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
		
		//we don't want to divide by zero.. add a little correction (causes big positive numbers for scaling subtractions but thats ok)
		double correctForZero = 0.0000000001;
		double irHigh = ir.high == 0 ? correctForZero : ir.high;
		double irModr = ir.moderate == 0 ? correctForZero : ir.moderate;
		double irLow = ir.low == 0 ? correctForZero : ir.low;
		double irModf = ir.modifier == 0 ? correctForZero : ir.modifier;
		
		//alright.. let's test if 'high' is our scaling reference:
		int highScaleModrDiff = (int)Math.round(-(nrOfModerate-(nrOfHigh*(irModr/irHigh))));
		int highScaleLowDiff = (int)Math.round(-(nrOfLow-(nrOfHigh*(irLow/irHigh))));
		int highScaleModfDiff = (int)Math.round(-(nrOfModifier-(nrOfHigh*(irModf/irHigh))));
		
		//no? then check if we should scale on 'moderate'
		int modrScaleHighDiff = (int)Math.round(-(nrOfHigh-(nrOfModerate*(irHigh/irModr))));
		int modrScaleLowDiff = (int)Math.round(-(nrOfLow-(nrOfModerate*(irLow/irModr))));
		int modrScaleModfDiff = (int)Math.round(-(nrOfModifier-(nrOfModerate*(irModf/irModr))));
		
		//no? then check if we should scale on 'low'
		int lowScaleHighDiff = (int)Math.round(-(nrOfHigh-(nrOfLow*(irHigh/irLow))));
		int lowScaleModrDiff = (int)Math.round(-(nrOfModerate-(nrOfLow*(irModr/irLow))));
		int lowScaleModfDiff = (int)Math.round(-(nrOfModifier-(nrOfLow*(irModf/irLow))));
		
		//no? then check if we should scale on 'modifier'
		int modfScaleHighDiff = (int)Math.round(-(nrOfHigh-(nrOfModifier*(irHigh/irModf))));
		int modfScaleModrDiff = (int)Math.round(-(nrOfModerate-(nrOfModifier*(irModr/irModf))));
		int modfScaleLowDiff = (int)Math.round(-(nrOfLow-(nrOfModifier*(irLow/irModf))));
		
//		System.out.println("scaling subtractions for HIGH: moderate=" + highScaleModrDiff + ", low" + highScaleLowDiff + ", modifier" + highScaleModfDiff);
//		System.out.println("scaling subtractions for MODERATE: high=" + modrScaleHighDiff + ", low=" + modrScaleLowDiff + ", modifier=" + modrScaleModfDiff);
//		System.out.println("scaling subtractions for LOW: high=" + lowScaleHighDiff + ", moderate=" + lowScaleModrDiff + ", modifier=" + lowScaleModfDiff);
//		System.out.println("scaling subtractions for MODIFIER: high=" + modfScaleHighDiff + ", moderate=" + modfScaleModrDiff + ", low=" + modfScaleLowDiff);
		
		int removeFromHigh = 0;
		int removeFromModerate = 0;
		int removeFromLow = 0;
		int removeFromModifier = 0;
		
		if(highScaleModrDiff < 0 && highScaleLowDiff < 0 && highScaleModfDiff < 0)
		{
//			System.out.println("we must scale on HIGH impact using " + highScaleModrDiff + ", " + highScaleLowDiff + ", " + highScaleModfDiff);
			removeFromModerate = -highScaleModrDiff;
			removeFromLow = -highScaleLowDiff;
			removeFromModifier = -highScaleModfDiff;
			
		}
		else if(modrScaleHighDiff < 0 && modrScaleLowDiff < 0 && modrScaleModfDiff < 0)
		{
//			System.out.println("we must scale on MODERATE impact using " + modrScaleHighDiff + ", " + modrScaleLowDiff + ", " + modrScaleModfDiff);
			removeFromHigh = -modrScaleHighDiff;
			removeFromLow = -modrScaleLowDiff;
			removeFromModifier = -modrScaleModfDiff;
		}
		else if(lowScaleHighDiff < 0 && lowScaleModrDiff < 0 && lowScaleModfDiff < 0)
		{
//			System.out.println("we must scale on LOW impact using " + lowScaleHighDiff + ", " + lowScaleModrDiff + ", " + lowScaleModfDiff);
			removeFromHigh = -lowScaleHighDiff;
			removeFromModerate = -lowScaleModrDiff;
			removeFromModifier = -lowScaleModfDiff;
		}
		else if(modfScaleHighDiff < 0 && modfScaleModrDiff < 0 && modfScaleLowDiff < 0)
		{
//			System.out.println("we must scale on MODIFIER impact using " + modfScaleHighDiff + ", " + modfScaleModrDiff + ", " + modfScaleLowDiff);
			removeFromHigh = -modfScaleHighDiff;
			removeFromModerate = -modfScaleModrDiff;
			removeFromLow = -modfScaleLowDiff;
		}
		else
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
	
}
