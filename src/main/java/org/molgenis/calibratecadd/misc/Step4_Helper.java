package org.molgenis.calibratecadd.misc;

import java.util.ArrayList;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import org.apache.commons.math3.stat.descriptive.rank.Median;
import org.molgenis.calibratecadd.structs.ImpactRatios;
import org.molgenis.calibratecadd.structs.VariantIntersectResult;
import org.molgenis.data.Entity;

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
							System.out.println("match on chrom/pos/ref/alt:");
							System.out.println(exacVariant.toString());
							System.out.println(clinvarVariant.toString());
							
							//edit the variant if we have >1 alt so we only keep this particular alt allele and matching AF field
							if(altSplit.length > 1)
							{
								exacVariant.set("ALT", alt);
								exacVariant.set("AF", exacVariant.getString("AF").split(",", -1)[altIndex]);
								System.out.println("CHANGED TO: ");
								System.out.println(exacVariant.toString());
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



	public static List<Entity> filterExACvariantsByMAF(List<Entity> inExACOnly, double medianMAF)
	{
		List<Entity> res = new ArrayList<Entity>();
		filterVariants:
		for(Entity exacVariant : inExACOnly)
		{
			String[] altSplit = exacVariant.getString("ALT").split(",", -1);
			for(int altIndex = 0; altIndex < altSplit.length; altIndex++)
			{
				String alt = altSplit[altIndex];
				double maf = Double.parseDouble(exacVariant.getString("AF").split(",",-1)[altIndex]);
				int ac = Integer.parseInt(exacVariant.getString("AC_Adj").split(",",-1)[altIndex]);
		
				boolean keep = false;
				
				//the clinvar variants were all 'singletons', so we only select singletons from exac
				if(maf == 0.0 && ac == 1)
				{
					keep = true;
				}
				else if(maf <= medianMAF)
				{
					keep = true;
				}
				
				if(keep)
				{
					//'edit' variant for this particular alt allele
					if(altSplit.length > 1)
					{
						exacVariant.set("ALT", alt);
						exacVariant.set("AF", maf);
						exacVariant.set("AC_Adj", ac);
					}
					res.add(exacVariant);
					
					//don't consider other alt alleles, just stick with the one we found first and continue
					continue filterVariants;
				}
			}
		}
		return res;
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
		
		int total = nrOfHigh + nrOfModerate + nrOfLow + nrOfModifier;
		
		ImpactRatios ir = new ImpactRatios(high, moderate, low, modifier)
		
		System.out.println("counts: high=" + nrOfHigh + ", modr=" + nrOfModerate + ", low=" + nrOfLow + ", modf=" + nrOfModifier);
		return null;
	}



	public static List<Entity> shapeExACvariantsByImpactRatios(List<Entity> exacFilteredByMAF, ImpactRatios ir)
	{
		// TODO Auto-generated method stub
		return null;
	}
}
