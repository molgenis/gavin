package org.molgenis.data.annotation.joeri282exomes.struct;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Set;

/**
 * Structure to store allele counts per gene, per patient group.
 * For each combination, we store counts for pathogenic and VOUS variants.
 * Then we store 3 values: total allele count, recessive genotype count, dominant genotype count.
 * 
 * So the format is:
 * 			group1			group2
 * gene1	2,0,2|23,5,13	5,1,3|42,10,22
 * gene2	etc
 * 
 * @author jvelde
 *
 */
public class GeneGroupsAlleleCountUtils
{
	public GeneGroupsAlleleCountUtils()
	{
		genesToGroupsToAlleleCounts = new HashMap<String, HashMap<String, AlleleCounts>>();
		groups = new HashSet<String>();
		
	}
	
	HashMap<String, HashMap<String, AlleleCounts>> genesToGroupsToAlleleCounts;
	private Set<String> groups;
	private PrintWriter pw;
	
	
	public void readFromFile()
	{
		
	}
	
	public void writeToFile(File writeTo) throws FileNotFoundException
	{
		pw = new PrintWriter(writeTo);
		
		for(String group : groups)
		{
			pw.print("\t" + group);
		}
		pw.print("\n");
		
		for(String gene : genesToGroupsToAlleleCounts.keySet())
		{
			StringBuffer printLine = new StringBuffer();
			printLine.append(gene);
			for(String group : groups)
			{
				if(genesToGroupsToAlleleCounts.get(gene).containsKey(group))
				{
					printLine.append("\t" + genesToGroupsToAlleleCounts.get(gene).get(group).toString());
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
	
	public void addToPathoAlleleTable(String gene, String group, int addAmount)
	{
		groups.add(group);
		addToTable(gene, group);
		AlleleCounts ac = genesToGroupsToAlleleCounts.get(gene).get(group);
		ac.pathoTotalAlleles += addAmount;
	}
	
	public void addToPathoRecessiveTable(String gene, String group, int addAmount)
	{
		addToTable(gene, group);
		AlleleCounts ac = genesToGroupsToAlleleCounts.get(gene).get(group);
		ac.pathoRecessiveGenotypes += addAmount;
	}
	
	public void addToPathoDominantTable(String gene, String group, int addAmount)
	{
		addToTable(gene, group);
		AlleleCounts ac = genesToGroupsToAlleleCounts.get(gene).get(group);
		ac.pathoDominantGenotypes += addAmount;
	}
	
	public void addToVOUSAlleleTable(String gene, String group, int addAmount)
	{
		addToTable(gene, group);
		AlleleCounts ac = genesToGroupsToAlleleCounts.get(gene).get(group);
		ac.vousTotalAlleles += addAmount;
	}
	
	public void addToVOUSRecessiveTable(String gene, String group, int addAmount)
	{
		addToTable(gene, group);
		AlleleCounts ac = genesToGroupsToAlleleCounts.get(gene).get(group);
		ac.vousRecessiveGenotypes += addAmount;
	}
	
	public void addToVOUSDominantTable(String gene, String group, int addAmount)
	{
		addToTable(gene, group);
		AlleleCounts ac = genesToGroupsToAlleleCounts.get(gene).get(group);
		ac.vousDominantGenotypes += addAmount;
	}
	
	private void addToTable(String gene, String group)
	{
		if(!genesToGroupsToAlleleCounts.containsKey(gene))
		{
			HashMap<String, AlleleCounts> newGene = new HashMap<String, AlleleCounts>();
			genesToGroupsToAlleleCounts.put(gene, newGene);
		}
		
		if(!genesToGroupsToAlleleCounts.get(gene).containsKey(group))
		{
			AlleleCounts newGroup = new AlleleCounts();
			genesToGroupsToAlleleCounts.get(gene).put(group, newGroup);
		}
	}
	
	
}
