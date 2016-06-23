package org.molgenis.cgd;

import java.io.*;
import java.util.*;
import java.util.zip.GZIPInputStream;

import org.molgenis.cgd.CGDEntry.generalizedInheritance;

public class LoadCGD {

	/**

	 X-linked dominant diseases & genes
	 Main source (after which followup searches to OMIM etc):

	 https://en.wikipedia.org/wiki/X-linked_dominant_inheritance#List_of_dominant_X-linked_diseases
	 https://en.wikipedia.org/wiki/Category:X-linked_dominant_disorders

	 PHEX	X-linked hypophosphatemia (Vitamin D resistant rickets)
	 FMR1	Fragile X Syndrome
	 MECP2	Rett syndrome
	 COL4A3/COL4A4	Most cases of Alport syndrome
	 IKBKG	Incontinentia pigmenti
	 -	Giuffrè–Tsukahara syndrome
	 PORCN	Goltz syndrome (Focal dermal hypoplasia)
	 - João-Gonçalves syndrome
	 -	Aicardi syndrome
	 ALAS2	X-linked dominant protoporphyria
	 -	Bazex–Dupré–Christol syndrome
	 NSDHL	CHILD syndrome
	 EFNB1	Craniofrontonasal dysplasia
	 MED12	Lujan–Fryns syndrome
	 BCOR	Oculofaciocardiodental syndrome

	 */

	public static List<String> xlinkedDominantGenes = Arrays.asList("PHEX", "FMR1", "MECP2", "COL4A3","COL4A4", "IKBKG", "PORCN", "ALAS2", "NSDHL", "EFNB1", "MED12", "BCOR");
	
	public static Map<String, CGDEntry> loadCGD(File cgdFile) throws IOException
	{

		InputStream fileStream = new FileInputStream(cgdFile);
		InputStream gzipStream = new GZIPInputStream(fileStream);
		Reader decoder = new InputStreamReader(gzipStream, "UTF-8");
		BufferedReader buffered = new BufferedReader(decoder);

		//HashMap<String, CGDEntry> cgd = new HashMap<String, CGDEntry>();
		Map<String, CGDEntry> cgd =  new TreeMap<String, CGDEntry>(String.CASE_INSENSITIVE_ORDER);

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
				inherMode = generalizedInheritance.DOMINANT_OR_RECESSIVE;
			}
			else if(split[4].contains("AR"))
			{
				inherMode = generalizedInheritance.RECESSIVE;
			}
			else if(split[4].contains("AD"))
			{
				inherMode = generalizedInheritance.DOMINANT;
			}
			else if(split[4].contains("BG"))
			{
				inherMode = generalizedInheritance.BLOODGROUP;
			}
			else if(split[4].contains("XL"))
			{
				if(xlinkedDominantGenes.contains(split[0]))
				{
					inherMode = generalizedInheritance.XL_DOMINANT;
				}
				else
				{
					inherMode = generalizedInheritance.XL_RECESSIVE;
				}
			}

			entry.setGeneralizedInheritance(inherMode);

			cgd.put(split[0], entry);

		});

		buffered.close();
        
        return cgd;
	}
	
	public static void main(String[] args) throws IOException {
		// TODO Auto-generated method stub

		Map<String, CGDEntry> cgd = LoadCGD.loadCGD(new File("/Users/joeri/github/gavin/data/other/CGD_1jun2016.txt.gz"));

		for( String key : cgd.keySet())
		{
			System.out.println(cgd.get(key).getGene() + " - " + cgd.get(key).getManifestationCategories() + " - " + cgd.get(key).getManifestationCategoriesList().toString());
		}

	}

}
