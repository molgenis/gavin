package org.molgenis.calibratecadd;

import java.io.File;
import java.io.IOException;
import java.util.HashMap;
import java.util.List;

import org.molgenis.data.annotator.tabix.TabixReader;
import org.molgenis.data.annotator.tabix.TabixReader.Iterator;

public class GetMatchingVariantsFromExAC
{

	// download @ ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz
	public static void main(String[] args) throws Exception
	{
		HashMap<String, List<ClinVarVariant>> geneToCVV = Step1_GetClinVarPathogenic.getAsMap(new File(args[0]));
		GetMatchingVariantsFromExAC.getAsMap(new File(args[1]), geneToCVV);
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
