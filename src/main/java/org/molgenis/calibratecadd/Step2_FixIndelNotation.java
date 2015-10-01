package org.molgenis.calibratecadd;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.HashMap;
import java.util.List;
import java.util.Scanner;

import org.molgenis.data.Entity;
import org.molgenis.data.annotator.tabix.TabixVcfRepository;

public class Step2_FixIndelNotation
{

	/**
	 * Uses:
	 * [0] file produced in step 1
	 * [1] ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/clinvar.vcf.gz
	 * [2] output file
	 * 
	 * try to fix 'na' and '-'
	 * 
	 * note: http://www.ncbi.nlm.nih.gov/clinvar/docs/ftp_primer/
	 * "At present, ClinVar's VCF file is limited to records that have been assigned rs#."
	 * 
	 * @param args
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException
	{
		Scanner s = new Scanner(new File(args[0]));
		PrintWriter pw = new PrintWriter(new File(args[2]));
		
		TabixVcfRepository clinvarVcf = new TabixVcfRepository(new File(args[1]), "clinvar");
		
		int fixes = 0;
		int failedFixes = 0;
		
		//failed fixed due to missing the target: 627 @ 10, 605 @ 20, 604 @ 25+
		int fixWindowSearchSize = 25;
		
		String line = null;
		
		HashMap<String, String> posToLine = new HashMap<String, String>();
		
		while(s.hasNextLine())
		{
			line = s.nextLine();
			if(line.startsWith("#"))
			{
				continue;
			}
			
			String[] split = line.split("\t", -1);
			
			if(split[2].equals("-1")) {
				//no RS id, so not fixable (we don't expect it in the VCF file)
				continue;
			}
			
			//deletions will have '-' as alt, insertions '-' as ref
			//sometimes also 'na' is used, e.g. "na	T" for a delCTinsA variant
			//so we check them all, and if none are found, continue
			if(!(split[3].equals("na") || split[4].equals("na") || split[3].equals("-") || split[4].equals("-")))
			{
				continue;
			}
			
//			System.out.println("attempting to fix: " + line);
			
			
			
			
			//NOTE: deletions will be found some position before!
			//e.g. "47705505 T -" will be in the VCF as "47705504 AT A" !
			// "47635557 AG -" will be at "47635554	rs63749848 CAG C"
			//tricky? yes... but we can grab a window and match by RS id :-)
			
			List<Entity> records = clinvarVcf.query(split[0], Long.parseLong(split[1])-fixWindowSearchSize, Long.parseLong(split[1])+fixWindowSearchSize);
			
			
			boolean match = false;
			for(Entity e : records)
			{
				String id = e.getString("ID").replace("rs", "");
				String id2 = split[2];
			//	System.out.println(id + " vs " + id2);
				if(id.equals(id2))
				{
					match = true;
					fixes++;
//					System.out.println("match! " + id);
				}
			}
			
			if(!match)
			{
				failedFixes++;
//				System.out.println("NO MATCH FOR " + line);
			}
			
			
			//now, get it from the ClinVar-VCF and match by RS id
			//there can be multiple variants on a position, e.g. 11:108205721, but using RS we can still get the correct one!
			
//			
//	//		if(!split[2].equals("-1")) { System.out.println("!!!!! " + line); }
//			String posKey = split[0] + "_" + split[1] + "_" + split[3] + "_" + split[4];
//			
//			if(posToLine.containsKey(posKey))
//			{
//		//		System.out.println("ALREADY CONTAINS " + posKey + " --> " + line);
//			}
//			posToLine.put(posKey, line);
//			
			
			
		}
		
		System.out.println("SUCCEEDED FIXES: " + fixes);
		System.out.println("FAILED FIXES: " + failedFixes);

	}

}
