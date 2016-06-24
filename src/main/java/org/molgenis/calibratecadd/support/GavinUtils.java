package org.molgenis.calibratecadd.support;

import java.io.File;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Scanner;
import java.util.Set;

import org.molgenis.data.Entity;

import org.molgenis.data.annotation.entity.impl.gavin.GavinEntry;
import org.molgenis.data.annotation.entity.impl.snpEff.Impact;


public class GavinUtils
{
	HashMap<String, GavinEntry> geneToEntry = new HashMap<String, GavinEntry>();
	
	public GavinUtils(File ccgg) throws Exception
	{	
		Scanner s = new Scanner(ccgg);
		
		//skip header
		s.nextLine();
		
		String line;
		while(s.hasNextLine())
		{
			line = s.nextLine();

			GavinEntry e = new GavinEntry(line);
			geneToEntry.put(e.gene, e);
		}
		
	}

	public HashMap<String, GavinEntry> getGeneToEntry()
	{
		return geneToEntry;
	}

	public GavinEntry.Category getCategory(String gene)
	{
		return geneToEntry.get(gene).category;
	}
	
	public boolean contains(String gene)
	{
		return geneToEntry.containsKey(gene) ? true : false;
	}
	
	public static Set<String> getGenesFromAnn(String ann) throws Exception
	{
		if(ann == null){
			return null;
		}
		Set<String> genes = new HashSet<String>();
		String[] annSplit = ann.split(",", -1);
		for(String oneAnn : annSplit)
		{
			String[] fields = oneAnn.split("\\|", -1);
			String gene = fields[3];
			genes.add(gene);
		}
		if(genes.size() == 0)
		{
			throw new Exception("No genes for " + ann);
		}
		return genes;
	}
	
	public static Double getInfoForAllele(Entity record, String infoField, String altAllele) throws Exception
	{
		String info_STR = record.get(infoField) == null ? null : record.get(infoField).toString();
		if(info_STR == null)
		{
			return null;
		}
		String[] alts = record.getString("ALT").split(",", -1);
		String[] info_split = info_STR.split(",", -1);
	
		if(alts.length != info_split.length)
		{
			throw new Exception("length of alts not equal to length of info field for " + record);
		}
		
		for (int i = 0; i < alts.length; i++)
		{
			if(alts[i].equals(altAllele))
			{
				return  (info_split[i] != null && !info_split[i].equals(".")) ? Double.parseDouble(info_split[i]) : null;
			}
		}
		return null;
	}
	
	public static Impact getImpact(String ann, String gene, String allele) throws Exception
	{
		//get the right annotation entry that matches both gene and allele
		String findAnn = getAnn(ann, gene, allele);
		if(findAnn == null)
		{
			System.out.println("WARNING: failed to get impact for gene '"+gene+"', allele '"+allele+"' in " + ann);
			return null;
		}
		else
		{
			//from the right one, get the impact
			String[] fields = findAnn.split("\\|", -1);
			String impact = fields[2];
			return Impact.valueOf(impact);
		}
	}

	public static String getTranscript(String ann, String gene, String allele) throws Exception
	{
		//get the right annotation entry that matches both gene and allele
		String findAnn = getAnn(ann, gene, allele);
		if(findAnn == null)
		{
			System.out.println("WARNING: failed to get impact for gene '"+gene+"', allele '"+allele+"' in " + ann);
			return null;
		}
		else
		{
			//from the right one, get the impact
			String[] fields = findAnn.split("\\|", -1);
			String transcript = fields[6];
			return transcript;
		}
	}
	
	public static String getAnn(String ann, String gene, String allele) throws Exception
	{
		String[] annSplit = ann.split(",", -1);
		for(String oneAnn : annSplit)
		{
			String[] fields = oneAnn.split("\\|", -1);
			String geneFromAnn = fields[3];
			if(!gene.equals(geneFromAnn))
			{
				continue;
			}
			String alleleFromAnn = fields[0];
			if(!allele.equals(alleleFromAnn))
			{
				continue;
			}
			return oneAnn;
		}
		System.out.println("WARNING: annotation could not be found for " + gene + ", allele=" + allele + ", ann=" + ann);
		return null;
	}

	public static void main(String[] args) throws Exception
	{
		File ccgg = new File(args[0]);
		new GavinUtils(ccgg);

	}

}
