package org.molgenis.data.annotation.joeri282exomes;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Scanner;
import java.util.Set;

import org.molgenis.data.Entity;
import org.molgenis.data.annotation.cmd.CommandLineAnnotatorConfig;
import org.molgenis.data.annotation.entity.impl.SnpEffAnnotator.Impact;
import org.molgenis.data.vcf.VcfRepository;
import org.springframework.context.annotation.AnnotationConfigApplicationContext;
import org.springframework.stereotype.Component;

@Component
public class CalibratedExomeAnalysis
{

	private HashMap<String, String> sampleToGroup;
//	private HashMap<String, List<CandidateVariant>> sampleToPathogenicVariants;
	
	private List<CandidateVariant> pathogenicVariants;
	
	public void go(File vcfFile, File patientGroups, File ccggFile) throws Exception
	{
		pathogenicVariants = new ArrayList<CandidateVariant>();
		
		loadSampleToGroup(patientGroups);

		CCGGUtils ccgg = new CCGGUtils(ccggFile);
		
		VcfRepository vcf = new VcfRepository(vcfFile, "vcf");
		Iterator<Entity> it = vcf.iterator();

		while (it.hasNext())
		{
			Entity record = it.next();

			String filter = record.getString("FILTER");

			if (!filter.equals("PASS"))
			{
				continue;
			}
			
		//	System.out.println(record.toString());
//			String chr = record.getString("#CHROM");
//			String pos = record.getString("POS");
//			String ref = record.getString("REF");
			String altStr = record.getString("ALT");
			String exac_af_STR = record.get("EXAC_AF") == null ? null : record.get("EXAC_AF").toString();
			String exac_ac_hom_STR = record.get("EXAC_AC_HOM") == null ? null : record.get("EXAC_AC_HOM").toString();
			String exac_ac_het_STR = record.get("EXAC_AC_HET") == null ? null : record.get("EXAC_AC_HET").toString();
			
			String ann = record.getString("ANN");

			String cadd_STR = record.get("CADD_SCALED") == null ? null : record.get("CADD_SCALED").toString();
			
			String[] alts = altStr.split(",", -1);

			
			String[] exac_af_split = new String[alts.length];
			if (exac_af_STR != null)
			{
				exac_af_split = exac_af_STR.split(",", -1);
			}
			String[] exac_ac_hom_split = new String[alts.length];
			if (exac_ac_hom_STR != null)
			{
				exac_ac_hom_split = exac_ac_hom_STR.split(",", -1);
			}

			String[] exac_ac_het_split = new String[alts.length];
			if (exac_ac_het_STR != null)
			{
				exac_ac_het_split = exac_ac_het_STR.split(",", -1);
			}
			String[] cadd_split = new String[alts.length];
			if (cadd_STR != null)
			{
				cadd_split = cadd_STR.split(",", -1);
			}
			
			Set<String> genes = getGenesFromAnn(ann);
			
			
			/**
			 * Iterate over alternatives, if applicable multi allelic example: 1:1148100-1148100
			 */
			for (int i = 0; i < alts.length; i++)
			{
				String alt = alts[i];
				
				double exac_af = (exac_af_split[i] != null && !exac_af_split[i].equals(".")) ? Double.parseDouble(exac_af_split[i]) : 0;
				Double cadd = (cadd_split[i] != null && !cadd_split[i].equals(".")) ? Double.parseDouble(cadd_split[i]) : null;
				int exac_hom = exac_ac_hom_split[i] != null && !exac_ac_hom_split[i].equals(".") ? Integer.parseInt(exac_ac_hom_split[i]) : 0;
				int exac_het = exac_ac_het_split[i] != null && !exac_ac_het_split[i].equals(".") ? Integer.parseInt(exac_ac_het_split[i]) : 0;
				
				for(String gene : genes)
				{
					if(!ccgg.geneToEntry.containsKey(gene))
					{
						continue;
					}
					Impact impact = getImpact(ann, gene, alt);

					
					Judgment j = ccgg.classifyVariant(gene, exac_af, impact, cadd);
					
					if(j.classification.equals(Judgment.Classification.Pathogenic))
					{
						List<String> samples = findInterestingSamples(record, i, exac_het, exac_hom);
						if(samples.size() > 0)
						{
							CandidateVariant cv = new CandidateVariant(record, alt, gene, j, samples);
							pathogenicVariants.add(cv);
							System.out.println("added patho variant: " + cv.toString());
						}
						
					}

				}
				
				
			}

		}
		vcf.close();
	}
	
	public List<String> findInterestingSamples(Entity record, int altIndex, int exacHets, int exacHoms)
	{
		List<String> homSampleIds = new ArrayList<String>();
		List<String> hetSampleIds = new ArrayList<String>();
		Iterable<Entity> sampleEntities = record.getEntities(VcfRepository.SAMPLES);
		//first alt has index 0, but we want to check 0/1, 1/1 etc. So +1.
		altIndex = altIndex + 1;
		
		int hetCount = 0;
		int homCount = 0;
		
		for (Entity sample : sampleEntities)
		{
			String genotype = sample.get("GT").toString();

			if (genotype.equals("./."))
			{
				continue;
			}

			// quality filter: we want depth X or more
			int depthOfCoverage = Integer.parseInt(sample.get("DP").toString());
			if (depthOfCoverage < 10)
			{
				continue;
			}

			
			if (genotype.equals("0/" + altIndex) || genotype.equals(altIndex + "/0")
					|| genotype.equals("0|" + altIndex) || genotype.equals(altIndex + "|0") )
			{
				//interesting!
				String sampleName = sample.get("ORIGINAL_NAME").toString();
				hetSampleIds.add(sampleName);
				hetCount++;
			}
			
			
			if ( genotype.equals(altIndex + "/" + altIndex) || genotype.equals(altIndex + "|" + altIndex) )
			{
				//interesting!
				String sampleName = sample.get("ORIGINAL_NAME").toString();
				homSampleIds.add(sampleName);
				homCount++;
			}

			
		}
		
		//TODO
		//do something based on exac genotype counts? e.g. delete heterozygous set when > 10 GTC or so
//		if(exacHoms > 10)
//		{
//			return new ArrayList<String>();
//		}
		
//		System.out.println("homs: " + homCount +" , hets: " + hetCount);
//		System.out.println("exachoms: " + exacHoms +" , exachets: " + exacHets);
		List<String> sampleIds = new ArrayList<String>();
		sampleIds.addAll(homSampleIds);
		sampleIds.addAll(hetSampleIds);
		return sampleIds;
	}
	
	public Impact getImpact(String ann, String gene, String allele) throws Exception
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
			String impact = fields[2];
			return Impact.valueOf(impact);
		}
		System.out.println("warning: impact could not be determined for " + gene + ", allele=" + allele + ", ann=" + ann);
		return null;
	}
	
	public Set<String> getGenesFromAnn(String ann)
	{
		Set<String> genes = new HashSet<String>();
		String[] annSplit = ann.split(",", -1);
		for(String oneAnn : annSplit)
		{
			String[] fields = oneAnn.split("\\|", -1);
//			String impact = fields[2];
//			String effect = fields[1];
//			String cDNA = fields[9];
//			String aaChange = fields[10];
			String gene = fields[3];
			genes.add(gene);
		}
		return genes;
	}

	public static void main(String[] args) throws Exception
	{
		// configureLogging();

		// See http://stackoverflow.com/questions/4787719/spring-console-application-configured-using-annotations
		AnnotationConfigApplicationContext ctx = new AnnotationConfigApplicationContext("org.molgenis.data.annotation");
		ctx.register(CommandLineAnnotatorConfig.class);
		CalibratedExomeAnalysis main = ctx.getBean(CalibratedExomeAnalysis.class);
		main.run(args);
		ctx.close();
	}
	

	private void loadSampleToGroup(File patientGroups) throws FileNotFoundException
	{
		HashMap<String, String> sampleToGroup = new HashMap<String, String>();
		Scanner s = new Scanner(patientGroups);
		String line = null;
		while (s.hasNextLine())
		{
			line = s.nextLine();
			String[] lineSplit = line.split("\t", -1);
			sampleToGroup.put(lineSplit[0], lineSplit[1]);
		}
		s.close();
		this.sampleToGroup = sampleToGroup;
	}

	public void run(String[] args) throws Exception
	{
		//TODO: first pass, print list of variants that miss a CADD score
		//we can do SNVs easy but not indels and other 'non predictable' changes
		//send off to CADD webservice, and give the output so we use this in second pass
		
		if (!(args.length == 3))
		{
			throw new Exception("Must supply 3 arguments");
		}

		File vcfFile = new File(args[0]);
		if (!vcfFile.isFile())
		{
			throw new Exception("Input VCF file does not exist or directory: " + vcfFile.getAbsolutePath());
		}

		File patientGroups = new File(args[1]);
		if (!patientGroups.isFile())
		{
			throw new Exception("Patient groups file does not exist or directory: " + patientGroups.getAbsolutePath());
		}
		
		File ccggFile = new File(args[2]);
		if (!ccggFile.isFile())
		{
			throw new Exception("CCGG file does not exist or directory: " + ccggFile.getAbsolutePath());
		}

		CalibratedExomeAnalysis cf = new CalibratedExomeAnalysis();
		cf.go(vcfFile, patientGroups, ccggFile);

	}

}
