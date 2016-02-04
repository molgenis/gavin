package org.molgenis.data.annotation.joeri282exomes.struct;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Scanner;
import java.util.Set;

/**
 * Structure to store allele count objects per gene, per patient group.
 * Has helper functions to easily add counts to each gene/patient combination.
 * Can also read and write from file.
 * 
 * @author jvelde
 *
 */
public class GeneGroupsAlleleCountUtils
{
	public GeneGroupsAlleleCountUtils()
	{
		this.genesToGroupsToAlleleCounts = new HashMap<String, HashMap<String, AlleleCounts>>();
		this.groups = new HashSet<String>();
		
	}
	
	private HashMap<String, HashMap<String, AlleleCounts>> genesToGroupsToAlleleCounts;
	private Set<String> groups;
	
	
	
	public HashMap<String, HashMap<String, AlleleCounts>> getGenesToGroupsToAlleleCounts()
	{
		return genesToGroupsToAlleleCounts;
	}

	public Set<String> getGroups()
	{
		return groups;
	}

	public void readFromFile(String fileLoc) throws Exception
	{
		File f = new File(fileLoc);
		Scanner s = new Scanner(f);
		String header = s.nextLine();
		String[] headerSplit = header.split("\t");
		for(String groupInHeader : headerSplit)
		{
			if(!groupInHeader.isEmpty())
			{
				this.groups.add(groupInHeader);
			}
		}
		
		String[] groups = this.getGroups().toArray(new String[this.getGroups().size()]);
		String line = null;
		while(s.hasNextLine())
		{
			line = s.nextLine();
			String[] lineSplit = line.split("\t");
			String gene = null;
			for(int i = 0; i < lineSplit.length; i++)
			{
				if(i == 0)
				{
					gene = lineSplit[i];
					continue;
				}
				
				String group = groups[i-1];
				String acString = lineSplit[i];
				
			//	System.out.println(gene + " - " + group + " - " + lineSplit[i]);
				
				AlleleCounts ac = new AlleleCounts(acString);
				
				
				addToTable(gene, group);
				
				this.genesToGroupsToAlleleCounts.get(gene).put(group, ac);
				
			}
		}
		
		
	}
	
	public void writeToFile(File writeTo) throws FileNotFoundException
	{
		PrintWriter pw = new PrintWriter(writeTo);
		
		for(String group : groups)
		{
			pw.print("\t" + group);
		}
		pw.print("\n");
		
		for(String gene : this.genesToGroupsToAlleleCounts.keySet())
		{
			StringBuffer printLine = new StringBuffer();
			printLine.append(gene);
			for(String group : groups)
			{
				if(this.genesToGroupsToAlleleCounts.get(gene).containsKey(group))
				{
					printLine.append("\t" + this.genesToGroupsToAlleleCounts.get(gene).get(group).toString());
				}
				else
				{
					printLine.append("\t" + AlleleCounts.printNull);
				}
				
			}
			pw.println(printLine.toString());
		}
		pw.flush();
		pw.close();
	}
	
	
	

	@Override
	public String toString()
	{
		return "GeneGroupsAlleleCountUtils [genesToGroupsToAlleleCounts=" + genesToGroupsToAlleleCounts + ", groups="
				+ groups + "]";
	}

	public void addToActingDominant(String gene, String group, int addAmount)
	{
		addToTable(gene, group);
		AlleleCounts ac = this.genesToGroupsToAlleleCounts.get(gene).get(group);
		ac.setActingDominant(ac.getActingDominant() + addAmount);
		this.genesToGroupsToAlleleCounts.get(gene).put(group, ac);
	}
	
	public void addToActingAdditive(String gene, String group, int addAmount)
	{
		
		addToTable(gene, group);
		AlleleCounts ac = this.genesToGroupsToAlleleCounts.get(gene).get(group);
		ac.setActingAdditive(ac.getActingAdditive() + addAmount);
	}
	
	public void addToActingRecessive(String gene, String group, int addAmount)
	{
		addToTable(gene, group);
		AlleleCounts ac = this.genesToGroupsToAlleleCounts.get(gene).get(group);
		ac.setActingRecessive(ac.getActingRecessive() + addAmount);
	}
	
	public void addToNonActingDominant(String gene, String group, int addAmount)
	{
		addToTable(gene, group);
		AlleleCounts ac = this.genesToGroupsToAlleleCounts.get(gene).get(group);
		ac.setNonActingDominant(ac.getNonActingDominant() + addAmount);
	}
	
	public void addToNonActingAdditive(String gene, String group, int addAmount)
	{
		addToTable(gene, group);
		AlleleCounts ac = this.genesToGroupsToAlleleCounts.get(gene).get(group);
		ac.setNonActingAdditive(ac.getNonActingAdditive() + addAmount);
	}
	
	public void addToNonActingRecessive(String gene, String group, int addAmount)
	{
		addToTable(gene, group);
		AlleleCounts ac = this.genesToGroupsToAlleleCounts.get(gene).get(group);
		ac.setNonActingRecessive(ac.getNonActingRecessive() + addAmount);
	}
	
	private void addToTable(String gene, String group)
	{
		groups.add(group);
		
		if(!this.genesToGroupsToAlleleCounts.containsKey(gene))
		{
			HashMap<String, AlleleCounts> newGene = new HashMap<String, AlleleCounts>();
			this.genesToGroupsToAlleleCounts.put(gene, newGene);
		}
		
		if(!this.genesToGroupsToAlleleCounts.get(gene).containsKey(group))
		{
			AlleleCounts newGroup = new AlleleCounts();
			this.genesToGroupsToAlleleCounts.get(gene).put(group, newGroup);
		}
		
		
	//	System.out.println("current AC for " + gene + " in group " + group + " = " + this.genesToGroupsToAlleleCounts.get(gene).get(group).toString());
	}
	
	
}
