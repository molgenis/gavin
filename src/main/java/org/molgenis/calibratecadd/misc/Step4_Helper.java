package org.molgenis.calibratecadd.misc;

import java.util.ArrayList;
import java.util.List;

import org.molgenis.calibratecadd.structs.ImpactRatios;
import org.molgenis.calibratecadd.structs.VariantIntersectResult;
import org.molgenis.data.Entity;

public class Step4_Helper
{

	public static VariantIntersectResult intersectVariants(List<Entity> exac, List<Entity> clinvar)
	{
		List<Entity> res = new ArrayList<Entity>();
		for (Entity exacVariant : exac)
		{
			for (Entity clinvarVariant : clinvar)
			{
				if (exacVariant.getString("#CHROM").equals(clinvarVariant.getString("#CHROM"))
						&& exacVariant.getString("POS").equals(clinvarVariant.getString("POS")))
				{

					System.out.println("match on chrom/pos:");
					System.out.println(exacVariant.toString());
					System.out.println(clinvarVariant.toString());

					if (exacVariant.getString("REF").equals(clinvarVariant.getString("REF"))
							&& exacVariant.getString("ALT").equals(clinvarVariant.getString("ALT")))
					{
						System.out.println("\talso match on ref/alt");
						System.out.println("\t" + exacVariant.toString());
						System.out.println("\t" + clinvarVariant.toString());
						// TODO match multiple alts!!!
						// MapEntity{#CHROM='2', ALT='T,A', POS='27693963', REF='C', FILTER='PASS', QUAL='7817.35',
						// ID='', INTERNAL_ID='2_27693963_C_T,A', INFO='null', SAMPLES_ENTITIES='[]'}
						// MapEntity{#CHROM='2', ALT='T', POS='27693963', REF='C', FILTER='null', QUAL='null',
						// ID='370540673', INTERNAL_ID='2_27693963_C_T', INFO='null'}

					}
				}
			}
		}
		
		
		//new VariantIntersectResult 
		
		return null;
	}
	
	
	
	public static double calculateMedianMAF(List<Entity> exacVariants)
	{
		return 0.0;
	}



	public static List<Entity> filterExACvariantsByMAF(List<Entity> inExACOnly, double medianMAF)
	{
		// TODO Auto-generated method stub
		return null;
	}



	public static ImpactRatios calculateImpactRatios(List<Entity> inClinVarOnly, List<Entity> inBoth)
	{
		// TODO Auto-generated method stub
		return null;
	}



	public static List<Entity> shapeExACvariantsByImpactRatios(List<Entity> exacFilteredByMAF, ImpactRatios ir)
	{
		// TODO Auto-generated method stub
		return null;
	}
}
