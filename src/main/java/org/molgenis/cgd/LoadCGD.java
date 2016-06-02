package org.molgenis.cgd;

import java.io.*;
import java.util.HashMap;
import java.util.Scanner;
import java.util.zip.GZIPInputStream;

import org.molgenis.cgd.CGDEntry.generalizedInheritance;

public class LoadCGD {

	
	
	public static HashMap<String, CGDEntry> loadCGD(File cgdFile) throws IOException
	{

		InputStream fileStream = new FileInputStream(cgdFile);
		InputStream gzipStream = new GZIPInputStream(fileStream);
		Reader decoder = new InputStreamReader(gzipStream, "UTF-8");
		BufferedReader buffered = new BufferedReader(decoder);

		HashMap<String, CGDEntry> cgd = new HashMap<String, CGDEntry>();

		buffered.lines().forEach(line -> {

			String[] split = line.split("\t", -1);

			CGDEntry entry = new CGDEntry(split[0], split[1], split[2], split[3], split[4], split[5], split[6], split[7], split[8], split[9], split[10], split[11]);

			// How to match these correctly? dozens of different terms, though most are "AR", "AD", etc.
			// However there are many combinations and exceptions.
			// Does AR take prevalence of AD or the other way around? be restrictive or loose here?

			// TODO: correctly!!

			generalizedInheritance inherMode = generalizedInheritance.OTHER;
			if(split[4].contains("AD") && split[4].contains("AR"))
			{
				inherMode = generalizedInheritance.DOM_OR_REC;
			}
			else if(split[4].contains("AR"))
			{
				inherMode = generalizedInheritance.RECESSIVE;
			}
			else if(split[4].contains("AD"))
			{
				inherMode = generalizedInheritance.DOMINANT;
			}
			else if(split[4].contains("XL"))
			{
				inherMode = generalizedInheritance.XLINKED;
			}

			entry.setGeneralizedInheritance(inherMode);

			cgd.put(split[0], entry);

		});

		buffered.close();
        
        return cgd;
	}
	
	public static void main(String[] args) throws IOException {
		// TODO Auto-generated method stub

		HashMap<String, CGDEntry> cgd = LoadCGD.loadCGD(new File("/Users/joeri/github/gavin/data/other/CGD_1jun2016.txt.gz"));

		for( String key : cgd.keySet())
		{
			System.out.println(cgd.get(key).getGene() + " - " + cgd.get(key).getManifestationCategories() + " - " + cgd.get(key).getManifestationCategoriesList().toString());
		}

	}

}
