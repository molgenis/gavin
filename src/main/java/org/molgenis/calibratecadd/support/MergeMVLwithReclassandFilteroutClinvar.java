package org.molgenis.calibratecadd.support;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

import org.molgenis.data.Entity;
import org.molgenis.data.MolgenisDataException;
import org.molgenis.data.vcf.VcfRepository;
import org.molgenis.data.vcf.utils.VcfUtils;
import org.molgenis.data.vcf.utils.VcfWriterUtils;

public class MergeMVLwithReclassandFilteroutClinvar
{
	//chrom_pos_ref_alt as key
	HashMap<String, Entity> mvlRecords = new HashMap<String, Entity>();
	
	//chrom_pos_ref_alt as key
	HashMap<String, Entity> reclassfRecords = new HashMap<String, Entity>();
	
	//chrom_pos_ref_alt
	List<String> clinvarPositions = new ArrayList<String>();
	
	public MergeMVLwithReclassandFilteroutClinvar(File mvlFile, File reclassfFile, File clinvarFile, File output) throws IOException
	{
		readFiles(mvlFile, reclassfFile, clinvarFile);
		

		System.out.println("MVL size: " + mvlRecords.size());
		System.out.println("Reclsf size: " + reclassfRecords.size());
		System.out.println("ClinVar size: " + clinvarPositions.size());
		
		HashMap<String, Entity> mergedRecords = mergeMvlWithReclassf();
		System.out.println("MVL merged with Reclsf size: " + mergedRecords.size());
		
		HashMap<String, Entity> mergedFilteredRecords = removeClinVar(mergedRecords);
		System.out.println("After removing ClinVar variants size: " + mergedFilteredRecords.size());
		
		print(mergedFilteredRecords, output);
	}
	
	public void print(HashMap<String, Entity> printMe, File output) throws MolgenisDataException, IOException
	{
		BufferedWriter pw = new BufferedWriter(new PrintWriter(output));
		for(String key : printMe.keySet())
		{
			VcfWriterUtils.writeToVcf(printMe.get(key), pw);
			pw.write('\n');
		}
		pw.flush();
		pw.close();
	}
	
	public HashMap<String, Entity> removeClinVar(HashMap<String, Entity> mergedRecords)
	{
		HashMap<String, Entity> filteredRecords = new HashMap<String, Entity>();
		
		for(String key : mergedRecords.keySet())
		{
			if(!clinvarPositions.contains(key))
			{
				filteredRecords.put(key, mergedRecords.get(key));
			}
		}
		return filteredRecords;
	}
	
	public HashMap<String, Entity> mergeMvlWithReclassf()
	{
		HashMap<String, Entity> mergedRecords = new HashMap<String, Entity>();
		
		List<String> inBoth = new ArrayList<String>();
		for(String key : mvlRecords.keySet())
		{
			if(reclassfRecords.containsKey(key))
			{
				inBoth.add(key);
				
				if(!reclassfRecords.get(key).getString("CLSF").equals(mvlRecords.get(key).getString("CLSF")))
				{
					System.out.println("classification required updating for " + key + "! " + mvlRecords.get(key).getString("CLSF") + " in MVL, " + reclassfRecords.get(key).getString("CLSF") + " in reclassified list");
				}
				//it shouldn't matter if we add from reclassfRecords or mvlRecords, but reclassfRecords will have the latest classification!
				mergedRecords.put(key, reclassfRecords.get(key));
			}
			else
			{
				mergedRecords.put(key, mvlRecords.get(key));
			}
			
		}
		
		for(String key : reclassfRecords.keySet())
		{
			if(!inBoth.contains(key))
			{
				mergedRecords.put(key, reclassfRecords.get(key));
			}
		}
		
		return mergedRecords;
	}
	
	public void readFiles(File mvlFile, File reclassfFile, File clinvarFile) throws IOException
	{
		VcfRepository mvlVcfRepo = new VcfRepository(mvlFile, "mvl");
		java.util.Iterator<Entity> vcfRepoIter = mvlVcfRepo.iterator();
		while (vcfRepoIter.hasNext())
		{
			Entity record = vcfRepoIter.next();
			mvlRecords.put(record.getString("#CHROM") + "_" + record.getString("POS") + "_" + record.getString("REF") + "_" + record.getString("ALT"), record);
		}
		mvlVcfRepo.close();
		
		VcfRepository reclassfVcfRepo = new VcfRepository(reclassfFile, "reclassifications");
		vcfRepoIter = reclassfVcfRepo.iterator();
		while (vcfRepoIter.hasNext())
		{
			Entity record = vcfRepoIter.next();
			reclassfRecords.put(record.getString("#CHROM") + "_" + record.getString("POS") + "_" + record.getString("REF") + "_" + record.getString("ALT"), record);
		}
		reclassfVcfRepo.close();
		
		VcfRepository clinvarVcfRepo = new VcfRepository(clinvarFile, "clinvar");
		vcfRepoIter = clinvarVcfRepo.iterator();
		while (vcfRepoIter.hasNext())
		{
			Entity record = vcfRepoIter.next();
			clinvarPositions.add(record.getString("#CHROM") + "_" + record.getString("POS") + "_" + record.getString("REF") + "_" + record.getString("ALT"));
		}
	}
	

	public static void main(String[] args) throws IOException
	{
		File mvlFile = new File(args[0]);
		File reclassfFile = new File(args[1]);
		File clinvarFile = new File(args[2]);
		File outputFile = new File(args[3]);
		new MergeMVLwithReclassandFilteroutClinvar(mvlFile, reclassfFile, clinvarFile, outputFile);

	}

}
