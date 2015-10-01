package org.molgenis.calibratecadd;

import java.io.File;
import java.io.IOException;
import java.util.HashMap;
import java.util.List;

import org.molgenis.calibratecadd.structs.ClinVarVariant;
import org.molgenis.data.annotator.tabix.TabixReader;
import org.molgenis.data.annotator.tabix.TabixReader.Iterator;

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
		
		
		
	}

	
	private static void getAsMap(File exac, HashMap<String, List<ClinVarVariant>> geneToCVV) throws IOException
	{
		System.out.println("loading matching exac variants..");
		
		for(String gene : geneToCVV.keySet())
		{
			String chrom = null;
			long leftMostPos = -1;
			long rightMostPos = -1;
			
			for(ClinVarVariant cvv : geneToCVV.get(gene))
			{
				if(cvv.pos > rightMostPos)
				{
					rightMostPos = cvv.pos;
				}
				
				if(cvv.pos < leftMostPos || leftMostPos == -1)
				{
					leftMostPos = cvv.pos;
				}
				if(chrom == null)
				{
					chrom = cvv.chrom;
				}
			}
			
			
			TabixReader tr = new TabixReader(exac.getAbsolutePath());
			
			Iterator it = tr.query(chrom + ":" + leftMostPos + "-" + rightMostPos);
			
			int count = 0;
			String next = null;
			while (it != null && (next = it.next()) != null)
    		{
			//	System.out.println(next);
				count++;
    		}
			
			
			System.out.println(gene + "\t"+leftMostPos + "\t" + rightMostPos + "\t" + "has " + "\t" + count );
			
		}
		
	}
	
	
	
	
}
