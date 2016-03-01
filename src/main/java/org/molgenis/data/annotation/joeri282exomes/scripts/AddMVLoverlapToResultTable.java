package org.molgenis.data.annotation.joeri282exomes.scripts;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Scanner;

import org.molgenis.data.Entity;
import org.molgenis.data.vcf.VcfRepository;

public class AddMVLoverlapToResultTable
{

	public static void main(String[] args) throws IOException
	{
		/** load MVL (VCF format)
		 * /Users/jvelde/Desktop/clinvarcadd/dataset_UMCGMVLsReworkedFeb2016/provenance/All_MVLs_reclassified.fix.snpeff.exac.cadd.vcf
			#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
			1	154574443	ADAR:c.675G>A	C	T	.	.	CADD=0.555756;CADD_SCALED=7.844;MVL=DYS;CLSF=LB;ANN=T|synonymous_variant|LOW|ADAR|ADAR|transcript|NM_001111.4|protein_coding|2/15|c.675G>A|p.Pro225Pro|917/6671|675/3681|225/1226||;EXAC_AF=1.647E-5;EXAC_AC_HOM=0;EXAC_AC_HET=2;
			etc
		 */
		VcfRepository vcf = new VcfRepository(new File("/Users/jvelde/Desktop/clinvarcadd/dataset_UMCGMVLsReworkedFeb2016/provenance/All_MVLs_reclassified.fix.snpeff.exac.cadd.vcf"), "vcf");
		Iterator<Entity> it = vcf.iterator();
		HashMap<String, String> chrPosRefAltToMVLdata = new HashMap<String, String>();
		while (it.hasNext())
		{
			Entity record = it.next();
			String chr = record.getString("#CHROM");
			String pos = record.getString("POS");
			String ref = record.getString("REF");
			String alt = record.getString("ALT");
			String mvlClassfc = record.getString("CLSF");
			String mvlName = record.getString("MVL");
			if(alt.contains(","))
			{
				System.out.println("WARNING: multi alt for" + record.toString());
			}
			String key = chr + "_" + pos + "_" + ref + "_" + alt;
			chrPosRefAltToMVLdata.put(key, "\t" + mvlName + "\t" + mvlClassfc);
		}
	
		/** parse result table file:
		 * /Users/jvelde/Desktop/joeri282exomes/joeri80dystonia/dystonia_allvariants_resulttable.txt
			Chrom	Pos	Ref	Alt	AffectingAlt	AffectedGene	cDNA	aaChange	Effect	ID	GoNL-AF	ExAC-AF	1000G-AF	CADD-score	Classification	Method	Reason
			4	17506013	G	A	A	QDPR	c.284C>T	p.Ala95Val	missense_variant		0	0	0	26.1	Pathogn	genomewide	Variant MAF of 0.0 is rare enough to be potentially pathogenic and the CADDscore of 26.1 is greater than a global threshold of 25.
			9	80022280	T	C	C	VPS13A	c.9400-169T>C		intron_variant	rs41289983	0.0963855421686747	0	0.0497204	4.006	Benign	calibrated	Variant is of 'modifier' impact, and therefore unlikely to be pathogenic. However, the variant MAF of 0.0 is less than the VPS13A pathogenic 95th percentile MAF of 1.5656E-5.
			etc
			
			Write out table with extra colums:
		 	[original header] \t MVL-Panel \t MVL-Classification
		 	[original line] \t DYS \t LB
		 */
		int overlap = 0;
		Scanner s2 = new Scanner(new File("/Users/jvelde/Desktop/joeri282exomes/joeri80dystonia/dystonia_allvariants_resulttable.txt"));
		String header = s2.nextLine();
		System.out.println(header + "\t" + "MVL-Panel" + "\t" + "MVL-Classification");
		while(s2.hasNextLine())
		{
			String line = s2.nextLine();
			String[] split = line.split("\t");
			String key = split[0] + "_" + split[1] + "_" + split[2] + "_" + split[4]; //4 is the acting allele out of multiple alleles
			if(chrPosRefAltToMVLdata.containsKey(key))
			{
				System.out.println(line + chrPosRefAltToMVLdata.get(key));
				overlap++;
			}
			else
			{
				System.out.println(line + "\t" + "-" + "\t" + "-");
			}
		}
		s2.close();

		System.out.println("\n\n###\noverlap: " + overlap + "\n###");
	}

}
