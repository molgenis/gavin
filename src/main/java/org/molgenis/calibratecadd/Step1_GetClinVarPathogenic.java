package org.molgenis.calibratecadd;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Scanner;

public class Step1_GetClinVarPathogenic
{


	// download @ ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz
	public static void main(String[] args) throws Exception
	{
		HashMap<String, List<ClinVarVariant>> cvv = Step1_GetClinVarPathogenic.getAsMap(new File(args[0]));
		PrintWriter pw = new PrintWriter(new File(args[1]));
		
		for(String chr : new String[]{"1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y", "MT"})
		{
			pw.println("##contig=<ID="+chr+">");
		}
		
		pw.println("##fileformat=VCFv4.1");
		pw.println("#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO");
		
		for(String gene : cvv.keySet())
		{
			for(ClinVarVariant v : cvv.get(gene))
			{
				pw.println(v.chrom + "\t" + v.pos + "\t" + v.id + "\t" + v.ref + "\t" + v.alt + "\t" + "." + "\t" + "." + "\t" + "CLINVAR=" + v.clinvarInfoToString());
			}
		}
		pw.flush();
		pw.close();
		
	}
	
	
	public static HashMap<String, List<ClinVarVariant>> getAsMap(File clinvarFile) throws FileNotFoundException
	{
		System.out.println("loading clinvar..");
		HashMap<String, List<ClinVarVariant>> res = new HashMap<String, List<ClinVarVariant>>();

		Scanner s = new Scanner(clinvarFile);

		String line = null;

		int lost = 0;
		int totalvariants = 0;
		
		while (s.hasNextLine())
		{
			line = s.nextLine();
			String[] lineSplit = line.split("\t", -1);

			// needs to be GRCh37
			String genomeBuild = lineSplit[12];
			if (!genomeBuild.equals("GRCh37"))
			{
				continue;
			}

			// needs to contain 'pathogenic'
			String clinsig = lineSplit[5];
			if (!clinsig.toLowerCase().contains("pathogenic"))
			{
				continue;
			}

			// eg. NM_005343.2(HRAS):c.37G>C (p.Gly13Arg)
			String name = lineSplit[2];
			// eg. HRAS
			String gene = lineSplit[4];

			// System.out.println(line);

			String geneFromName = null;
			if (name.contains("("))
			{
				geneFromName = name.substring(name.indexOf('(') + 1, name.indexOf(')'));
			}

			if (gene.equals("-") || gene.equals(""))
			{
				if(gene.equals("")) { System.out.println("override: '"+gene+"' to '"+geneFromName+"'"); }
				gene = geneFromName;
			}

			if (gene == null)
			{
				lost++;
				continue;
			}

			String chrom = lineSplit[13];
			String pos = lineSplit[14];
			String id = lineSplit[6];
			String ref = lineSplit[25];
			String alt = lineSplit[26];

			ClinVarVariant cvv = new ClinVarVariant(chrom, pos, id, ref, alt, name, gene, clinsig);
			
			if(res.containsKey(gene))
			{
				res.get(gene).add(cvv);
			}
			else
			{
				List<ClinVarVariant> cvvList = new ArrayList<ClinVarVariant>();
				cvvList.add(cvv);
				res.put(gene, cvvList);
			}
			
			totalvariants++;

		}
		
		s.close();

		System.out.println("..done, put " + totalvariants + " 'pathogenic' variants in " + res.size() + " genes, lost " + lost + " due to non-recoverable gene symbols");
		
		return res;

	}

}
