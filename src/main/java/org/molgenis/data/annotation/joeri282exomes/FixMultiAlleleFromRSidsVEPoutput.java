package org.molgenis.data.annotation.joeri282exomes;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.util.HashMap;
import java.util.Map;
import java.util.Scanner;

public class FixMultiAlleleFromRSidsVEPoutput
{

	File variBenchOriginalFile;
	File variBenchConvertedToVCF;
	PrintWriter pw;
	
	/**
	 * after using Ensembl VEP webservice to retrieve variants using RS id, you often get multiple alleles for each
	 * variant. while in the original file (e.g. VariBench PON-P2 neutral training/test sets) there is some information
	 * on which allele was the informative/interpreted one. the second file is custom, so meh, but this fix can be
	 * important as it concerns ~10% of the variants..
	 * @throws FileNotFoundException 
	 * 
	 **/
	public FixMultiAlleleFromRSidsVEPoutput(File variBenchOriginalFile, File variBenchConvertedToVCF, File outputFile) throws FileNotFoundException
	{
		if (!variBenchOriginalFile.isFile())
		{
			throw new FileNotFoundException("variBenchOriginalFile " + variBenchOriginalFile.getAbsolutePath()
					+ " does not exist or is directory");
		}
		if (!variBenchConvertedToVCF.isFile())
		{
			throw new FileNotFoundException("variBenchConvertedToVCF " + variBenchConvertedToVCF.getAbsolutePath()
					+ " does not exist or is directory");
		}
		if (outputFile.isFile())
		{
			System.out.println("Warning: output file " + outputFile.getAbsolutePath()
					+ " already exists, overwriting content!");
		}
		this.variBenchOriginalFile = variBenchOriginalFile;
		this.variBenchConvertedToVCF = variBenchConvertedToVCF;
		this.pw = new PrintWriter(outputFile);
	}
	
	public String swapBase(String in) throws Exception
	{
		if(in.equals("A")){ return "T"; }
		else if(in.equals("T")){ return "A"; }
		else if(in.equals("G")){ return "C"; }
		else if(in.equals("C")){ return "G"; }
		else{ throw new Exception("Cannot handle non base " + in); }
	}
	
	public void annotate() throws Exception
	{
		/**
		 * Header:
		 * Rs_id	Contig_Acc_version ....
		 * 
		 * Example lines:
		 * 11260579	NT_004350.19	700898	700898	37_1	NM_002978.2	705	705	2	C	CCC	G	CGC	G>C	NP_002969.2	P	Pro	179	R	Arg	p.Arg180Pro	p.R180P	Cluster 393	NP_002969.2	VariO:0128 variation affecting DNA|VariO:0129 DNA variation type|VariO:0322 DNA variation classification|VariO:0135 DNA chain variation|VariO:0136 DNA substitution|VariO:0316 transversion	VariO:0297 variation affecting RNA|VariO:0306 RNA variation type|VariO:0328 RNA variation classification|VariO:0312 RNA substitution|VariO:0316 transversion	VariO:0002 variation affecting protein|VariO:0012 protein variation type|VariO:0325 protein variation classification|VariO:0021 amino acid substitution
		 * 61758325	NT_010799.15	7425797	7425797	37_1	NM_002981.1	259	259	2	G	AGA	T	ATA	T>G	NP_002972.1	R	Arg	62	I	Ile	p.Ile63Arg	p.I63R	Cluster 5995	NP_002972.1	VariO:0128 variation affecting DNA|VariO:0129 DNA variation type|VariO:0322 DNA variation classification|VariO:0135 DNA chain variation|VariO:0136 DNA substitution|VariO:0316 transversion	VariO:0297 variation affecting RNA|VariO:0306 RNA variation type|VariO:0328 RNA variation classification|VariO:0312 RNA substitution|VariO:0316 transversion	VariO:0002 variation affecting protein|VariO:0012 protein variation type|VariO:0325 protein variation classification|VariO:0021 amino acid substitution
		 * 75433456	NT_010799.15	7349851	7349851	37_1	NM_002986.2	171	171	1	G	GTG	C	CTG	C>G	NP_002977.1	V	Val	10	L	Leu	p.Leu11Val	p.L11V	Cluster 217	NP_005614.2	VariO:0128 variation affecting DNA|VariO:0129 DNA variation type|VariO:0322 DNA variation classification|VariO:0135 DNA chain variation|VariO:0136 DNA substitution|VariO:0316 transversion	VariO:0297 variation affecting RNA|VariO:0306 RNA variation type|VariO:0328 RNA variation classification|VariO:0312 RNA substitution|VariO:0316 transversion	VariO:0002 variation affecting protein|VariO:0012 protein variation type|VariO:0325 protein variation classification|VariO:0021 amino acid substitution
		 */
		Scanner variBenchOriginalFileScanner = new Scanner(variBenchOriginalFile);
		//skip header
		variBenchOriginalFileScanner.nextLine();
		Map<String,String> altAlleleForRS = new HashMap<String,String>();
		while(variBenchOriginalFileScanner.hasNextLine())
		{
			//e.g.  '2	47630249	GA	G	1.373305	10.52'
			String line = variBenchOriginalFileScanner.nextLine();
	        if(line.startsWith("#"))
	        {
	            continue;
	        }
	        String[] split = line.split("\t");
           
            //so: "rs11260579" to "C", unless empty..
	        if(!split[9].isEmpty())
	        {
	        	 altAlleleForRS.put("rs"+split[0], split[9]);
	             //System.out.println("put: " + "rs"+split[0] + " to " + split[9]);
	        }
		}
		variBenchOriginalFileScanner.close();
		
		System.out.println("altAlleleForRS size " + altAlleleForRS.size());

		/**
		 * Now fix the VCF. Example line:
		 * 
		 * 1	1222267	rs11260579	G	A,C,T
		 * 
		 * So here, we must reduce this to G -> C
		 */
		Scanner vcfOutputScanner = new Scanner(variBenchConvertedToVCF);
		while(vcfOutputScanner.hasNextLine())
		{
			String line = vcfOutputScanner.nextLine();
	        if(line.startsWith("#"))
	        {
	        	pw.println(line);
	            continue;
	        }
	        String[] split = line.split("\t", -1);
	        String rsID = split[2];
	        String[] vcfAltAlleles = split[4].split(",", -1);
	        if(altAlleleForRS.containsKey(rsID))
	        {
	        	boolean success = false;
	        	for(String vcfAlt : vcfAltAlleles)
	        	{
	        		
	        		String base = altAlleleForRS.get(rsID);
	        		String swappedBase = swapBase(base);
	        		
	        		boolean checkDefaultBase = false;
	        		boolean checkSwappedBase = false;
	        		
	        		if(split[7].contains("|1|") && split[7].contains("|-1|"))
	        		{
	        			checkDefaultBase = true;
	        			checkSwappedBase = true;
	        		}
	        		else if(split[7].contains("|1|") )
	        		{
	        			checkDefaultBase = true;
	        		}
	        		else if(split[7].contains("|-1|") )
	        		{
	        			checkSwappedBase = true;
	        		}
	        		
	        		if((checkDefaultBase && vcfAlt.equals(base)) || (checkSwappedBase && vcfAlt.equals(swappedBase)))
	        		{
	        			StringBuffer fixedLine = new StringBuffer();
	        			for(int i = 0; i < split.length; i++)
	        			{
	        				if(i == 4)
	        				{
	        					fixedLine.append(vcfAlt + "\t");
	        				}
	        				else
	        				{
	        					fixedLine.append(split[i] + "\t");
	        				}
	        			}
	        			fixedLine.deleteCharAt(fixedLine.length()-1);
	        			pw.println(fixedLine.toString());
	        			success = true;
	        			break;
	        		}
	        	}
	        	if(!success)
	        	{
	        		System.out.println("WARNING: we have altAlleleForRS="+altAlleleForRS.get(rsID) + " for rsID="+rsID+", but no good match to VCF alt alleles in line: " + line);
	        	}
	        }
	        else
	        {
	        	if(vcfAltAlleles.length > 1)
	        	{
	        		//still ambiguous! don't print line
	        		System.out.println("WARNING: no alt info for rsID="+rsID+", but there are multiple alts in VCF line: " + line);
	        	}
	        	else {
	        		//no replacement and only 1 alt, just print this line
		        	pw.println(line);
	        	}
	        }
		}
		pw.flush();
	}
	
	public static void main(String[] args) throws Exception
	{
		File variBenchOriginalFile = new File(args[0]);
		File variBenchConvertedToVCF = new File(args[1]);
		File outputFile = new File(args[2]);
		FixMultiAlleleFromRSidsVEPoutput fix = new FixMultiAlleleFromRSidsVEPoutput(variBenchOriginalFile, variBenchConvertedToVCF, outputFile);
		fix.annotate();

	}

}
