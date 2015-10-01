package org.molgenis.calibratecadd;

import java.io.File;
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
	 * @throws Exception 
	 */
	public static void main(String[] args) throws Exception
	{
		Scanner s = new Scanner(new File(args[0]));
		PrintWriter pw = new PrintWriter(new File(args[2]));
		
		TabixVcfRepository clinvarVcf = new TabixVcfRepository(new File(args[1]), "clinvar");
		
		int fixes = 0;
		int failedFixes = 0;
		int lost = 0;
		int totalPassedVariants = 0;
		
		//failed fixed due to missing the target: 627 @ 10, 605 @ 20, 604 @ 25+
		int fixWindowSearchSize = 25;
		
		String line = null;
		
		HashMap<String, String> posToLine = new HashMap<String, String>();
		
		while(s.hasNextLine())
		{
			line = s.nextLine();
			if(line.startsWith("#"))
			{
				//header, just print and continue
				pw.println(line);
				continue;
			}
			
			String[] lineSplit = line.split("\t", -1);
			
			//deletions will have '-' as alt, insertions '-' as ref
			//sometimes also 'na' is used, e.g. "na	T" for a delCTinsA variant
			//so we check them all, and if none are found, assume its OK, print and continue
			//additional check: ref and alt the same! eg. rs10993994 (10:51549496-51549496) : T T in clinvar-tabdelim, T C in clinvar-vcf ...
			if(!(lineSplit[3].equals("na") || lineSplit[4].equals("na") || lineSplit[3].equals("-") || lineSplit[4].equals("-") || lineSplit[3].equals(lineSplit[4])))
			{
				totalPassedVariants++;
				pw.println(line);
				continue;
			}
			
			if(lineSplit[2].equals("-1")) {
				//no RS id, so not fixable (we don't expect it in the VCF file, and couldn't match it there either)
				// do not print!
				lost++;
				continue;
			}

			
			//now the fixing can begin!
			
			//NOTE: deletions will be found some position before!
			//e.g. "47705505 T -" will be in the VCF as "47705504 AT A" !
			// "47635557 AG -" will be at "47635554	rs63749848 CAG C"
			//tricky? yes... but we can grab a window and match by RS id :-)
			
			List<Entity> records = clinvarVcf.query(lineSplit[0], Long.parseLong(lineSplit[1])-fixWindowSearchSize, Long.parseLong(lineSplit[1])+fixWindowSearchSize);
			
			boolean match = false;
			for(Entity e : records)
			{
				String idFromClinVarVCF = e.getString("ID").replace("rs", "");
				String idFromOurClinVarPathoList = lineSplit[2];
				if(idFromClinVarVCF.equals(idFromOurClinVarPathoList))
				{
					if(match == true)
					{
						throw new Exception("double match?? " + line);
					}
					StringBuffer fixedLine = new StringBuffer();
					for(int i = 0; i < lineSplit.length; i++)
					{
						if(i == 3)
						{
							fixedLine.append(e.getString("REF") + "\t");
						}
						else if(i == 4)
						{
							fixedLine.append(e.getString("ALT") + "\t");
						}
						else
						{
							fixedLine.append(lineSplit[i] + "\t");
						}
//						if(e.getString("REF").equals(e.getString("ALT")))
//						{
//							System.out.println("WTF?? " + e.toString());
//						}
					}
					fixedLine.deleteCharAt(fixedLine.length()-1);
					pw.println(fixedLine.toString());
					match = true;
					fixes++;
					totalPassedVariants++;
				}
			}
			
			if(!match)
			{
				failedFixes++;
			}
			
		}
		
		pw.flush();
		pw.close();
		
		System.out.println("total variants seen:" + (totalPassedVariants+lost+failedFixes));
		System.out.println("total passed variants after checking & fixing: " + totalPassedVariants);
		System.out.println("total variants dropped (not fixable/fixed): " + (lost+failedFixes));
		System.out.println("variants dropped without attempting to fix (not possible): " + lost);
		System.out.println("total fixed attempted: " + (fixes+failedFixes));
		System.out.println("fixes succeeded: " + fixes);
		System.out.println("fixes failed: " + failedFixes);

	}

}
