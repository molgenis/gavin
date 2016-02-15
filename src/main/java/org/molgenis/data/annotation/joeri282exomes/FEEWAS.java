package org.molgenis.data.annotation.joeri282exomes;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Scanner;
import java.util.Set;

import org.molgenis.calibratecadd.support.LoadCADDWebserviceOutput;
import org.molgenis.calibratecadd.support.VariantClassificationException;
import org.molgenis.data.Entity;
import org.molgenis.data.annotation.cmd.CommandLineAnnotatorConfig;
import org.molgenis.data.annotation.entity.impl.SnpEffAnnotator.Impact;
import org.molgenis.data.annotation.joeri282exomes.struct.GeneGroupsAlleleCountUtils;
import org.molgenis.data.vcf.VcfRepository;
import org.springframework.context.annotation.AnnotationConfigApplicationContext;
import org.springframework.stereotype.Component;

/**
 * Fisher's Exact Exome-Wide Association Study
 * 
 * TODO
 * WORK IN PROGRESS PILOT
 * 
 * idea: classify variants, in context of allele & gene (so get rid of code for 'conflicts' and such)
 * then, count the number of genotypes/alleles in patient groups over genes
 * when we have gathered this data, use fisher exact test to find any genes with
 * overrepresented pathogenic (or VOUS) variants for specific patient groups
 * 
 * We partition the data per gene, so that the test is "1 group & other groups" vs. "nr of pathogenic genotypes & other genotypes".
 * This means the total amount of tests equals:
 * N genes * M groups * 3 inheritance modes * 2 types (patho/vous)
 * 
 * @author jvelde
 *
 */
@Component
public class FEEWAS
{

	private static final String MODE_ANALYSIS = "analysis";
	private static final String MODE_CREATECADDFILE = "create_file_for_cadd";
	
	private HashMap<String, String> sampleToGroup;
//	private HashMap<String, List<CandidateVariant>> sampleToPathogenicVariants;
	
	private List<CandidateVariant> pathogenicVariants;
	private List<Entity> noJudgmentVariants;
	private List<Entity> conflictingJudgmentVariants;
	private LinkedHashMap<String,String> geneLocs = new LinkedHashMap<String,String>();
	
	private int variantRefAltGenePatho = 0;
	private int variantRefAltGeneBenign = 0;
	private int variantSkippedBadFilter = 0;
	private int failedToClassify = 0;
	private int conflictingClassifications = 0;
	private int totalVariantSeen = 0;
	private int totalVariantRefAltGeneCombinationsSeen = 0;
	
	public void go(String mode, File vcfFile, File patientGroups, File ccggFile, File caddFile) throws Exception
	{
		long startTime = System.currentTimeMillis();
		pathogenicVariants = new ArrayList<CandidateVariant>();
		noJudgmentVariants = new ArrayList<Entity>();
		conflictingJudgmentVariants = new ArrayList<Entity>();
		
		loadSampleToGroup(patientGroups);
		GeneGroupsAlleleCountUtils geneGroupCountsPathogenic = new GeneGroupsAlleleCountUtils();
		GeneGroupsAlleleCountUtils geneGroupCountsVOUS = new GeneGroupsAlleleCountUtils();

		CCGGUtils ccgg = new CCGGUtils(ccggFile);
		
		VcfRepository vcf = new VcfRepository(vcfFile, "vcf");
		Iterator<Entity> it = vcf.iterator();
		
		//either print missing cadd scores to this file, or read from file to get them, depending on mode
		PrintWriter pw = null;
		if(mode.equals(MODE_CREATECADDFILE))
		{
			pw = new PrintWriter(caddFile);
		}
		HashMap<String, Double> caddScores = null;
		if(mode.equals(MODE_ANALYSIS))
		{
			caddScores = LoadCADDWebserviceOutput.load(caddFile);
		}
		
		while (it.hasNext())
		{
			if(totalVariantSeen % 1000 == 0)
			{
				System.out.println("Seen " + totalVariantSeen + " variants, in a total of " + totalVariantRefAltGeneCombinationsSeen + " ref/alt/gene combinations.");
			}
			totalVariantSeen++;
		
			
			Entity record = it.next();
			String filter = record.getString("FILTER");
			String chr = record.getString("#CHROM");
			String pos = record.getString("POS");
			String ref = record.getString("REF");
			String altStr = record.getString("ALT");
			String exac_af_STR = record.get("EXAC_AF") == null ? null : record.get("EXAC_AF").toString();
			String ann = record.getString("ANN");
			String cadd_STR = record.get("CADD_SCALED") == null ? null : record.get("CADD_SCALED").toString();
			String[] alts = altStr.split(",", -1);
			String[] exac_af_split = exac_af_STR != null ? exac_af_STR.split(",", -1) : null;
			String[] cadd_split = cadd_STR != null ? cadd_STR.split(",", -1) : null;
			Set<String> genes = CCGGUtils.getGenesFromAnn(ann);
			
			
			if (filter != null && !filter.equals("PASS"))
			{
				variantSkippedBadFilter++;
				continue;
			}
			
			/**
			 * Iterate over alternatives, if applicable multi allelic example: 1:1148100-1148100
			 */
			for (int i = 0; i < alts.length; i++)
			{
				String alt = alts[i];
				double exac_af = (exac_af_split != null && !exac_af_split[i].equals(".")) ? Double.parseDouble(exac_af_split[i]) : 0;
				
				/**
				 * Sort out CADD score for this variant/allele
				 */
				Double cadd = (cadd_split != null && !cadd_split[i].equals(".")) ? Double.parseDouble(cadd_split[i]) : null;
				if(cadd == null)
				{
					if(mode.equals(MODE_CREATECADDFILE))
					{
						String trimmedRefAlt = LoadCADDWebserviceOutput.trimRefAlt(ref, alt, "\t");
						pw.println(chr + "\t" + pos + "\t" + "." + "\t" + trimmedRefAlt);
					}
					else if(mode.equals(MODE_ANALYSIS))
					{
						String key = chr + "_" + pos + "_" + ref + "_" + alt;
						if(caddScores.containsKey(key))
						{
							cadd = caddScores.get(key);
						}
						else
						{
							String trimmedRefAlt = LoadCADDWebserviceOutput.trimRefAlt(ref, alt, "_");
							key = chr + "_" + pos + "_" + trimmedRefAlt;
							if(caddScores.containsKey(key))
							{
								cadd = caddScores.get(key);
							}
							else
							{
								System.out.println("WARNING: CADD score missing for " + chr + " " + pos + " " + ref + " " + alt + " ! (even when using trimmed key '"+key+"')");
							}
							
						}
					}
				}
				

				for(String gene : genes)
				{
					totalVariantRefAltGeneCombinationsSeen++;
					
					Impact impact = CCGGUtils.getImpact(ann, gene, alt);
					if(impact == null)
					{
						continue;
					}
					
					Judgment judgment = null;
					
					judgment = ccgg.classifyVariant(gene, exac_af, impact, cadd);
					
					if(judgment.classification.equals(Judgment.Classification.Pathogn))
					{
						
						countGroupSamples(record, i, geneGroupCountsPathogenic, gene);
						if(!geneLocs.containsKey(gene)) { geneLocs.put(gene, chr + "\t" + pos); }
				//		System.out.println("geneGroupCountsPathogenic=" + geneGroupCountsPathogenic.toString());
					}
					else if(judgment.classification.equals(Judgment.Classification.VOUS))
					{
						countGroupSamples(record, i, geneGroupCountsVOUS, gene);
						if(!geneLocs.containsKey(gene)) { geneLocs.put(gene, chr + "\t" + pos); }
					}

				
				}
				
				
			}

		}
		vcf.close();
		
		if(mode.equals(MODE_CREATECADDFILE))
		{
			pw.flush();
			pw.close();
		}
		geneGroupCountsPathogenic.setGeneLocs(geneLocs);
		geneGroupCountsVOUS.setGeneLocs(geneLocs);
		geneGroupCountsPathogenic.writeToFile(new File("/Users/jvelde/FEEWAS_Patho.tsv"));
		geneGroupCountsVOUS.writeToFile(new File("/Users/jvelde/FEEWAS_VOUS.tsv"));
		
		long endTime = System.currentTimeMillis();
		
		System.out.println("\n## COUNTS ##\n");
		System.out.println("total variant seen = " + totalVariantSeen);
		System.out.println("variants that did not pass QC and were skipped: " + variantSkippedBadFilter);
		System.out.println("total variants incl. altallele/gene combinations seen = " + totalVariantRefAltGeneCombinationsSeen);
		System.out.println("relevant combinations, classified benign = " + variantRefAltGeneBenign);
		System.out.println("relevant combinations, classified pathogenic = " + variantRefAltGenePatho);
		System.out.println("classification failed = " + failedToClassify);
		System.out.println("classification conflicting = " + conflictingClassifications);
		System.out.println("runtime: " + (endTime-startTime)/1000.0 + " seconds, " + totalVariantRefAltGeneCombinationsSeen/((endTime-startTime)/1000.0) + " interpretations per second");
	
	}
	
	public void countGroupSamples(Entity record, int altIndex, GeneGroupsAlleleCountUtils genegroupcounts, String gene)
	{
		Iterable<Entity> sampleEntities = record.getEntities(VcfRepository.SAMPLES);
		
		//first alt has index 0, but we want to check 0/1, 1/1 etc. So +1.
		altIndex = altIndex + 1;
		
		for (Entity sample : sampleEntities)
		{
			String genotype = sample.get("GT").toString();

			if (genotype.equals("./."))
			{
				continue;
			}

			if(sample.get("DP") == null)
			{
				continue;
			}
			
			// quality filter: we want depth X or more
			int depthOfCoverage = Integer.parseInt(sample.get("DP").toString());
			if (depthOfCoverage < 10)
			{
				continue;
			}

			String sampleName = sample.get("ORIGINAL_NAME").toString();
			String group = sampleToGroup != null ? sampleToGroup.get(sampleName) : "default";
			
		//	System.out.println("sample: " + sampleName + "\t" + genotype );
			
			/**
			 * Homozygous reference
			 */
			if (genotype.equals("0/0") || genotype.equals("0|0"))
			{
				genegroupcounts.addToNonActingDominant(gene, group, 1);
				genegroupcounts.addToNonActingAdditive(gene, group, 2);
				genegroupcounts.addToNonActingRecessive(gene, group, 1);
			}
			
			/**
			 * Heterozygous
			 */
			else if (genotype.equals("0/" + altIndex) || genotype.equals(altIndex + "/0")
					|| genotype.equals("0|" + altIndex) || genotype.equals(altIndex + "|0") )
			{
				genegroupcounts.addToNonActingAdditive(gene, group, 1);
				genegroupcounts.addToNonActingRecessive(gene, group, 1);
				genegroupcounts.addToActingDominant(gene, group, 1);
				genegroupcounts.addToActingAdditive(gene, group, 1);
				
			}
			
			/**
			 * Homozygous alternative
			 */
			else if ( genotype.equals(altIndex + "/" + altIndex) || genotype.equals(altIndex + "|" + altIndex) )
			{
				genegroupcounts.addToActingDominant(gene, group, 1);
				genegroupcounts.addToActingAdditive(gene, group, 2);
				genegroupcounts.addToActingRecessive(gene, group, 1);
			}

			else
			{
				//funky genotype, e.g. 1/2 or 2/3. found in hypervariable regions such as HLA a lot.
			}
		}

	}

	public static void main(String[] args) throws Exception
	{
		// configureLogging();

		// See http://stackoverflow.com/questions/4787719/spring-console-application-configured-using-annotations
		AnnotationConfigApplicationContext ctx = new AnnotationConfigApplicationContext("org.molgenis.data.annotation");
		ctx.register(CommandLineAnnotatorConfig.class);
		FEEWAS main = ctx.getBean(FEEWAS.class);
		main.run(args);
		ctx.close();
	}
	

	private void loadSampleToGroup(File patientGroups) throws FileNotFoundException
	{
		if(patientGroups.exists())
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
	}

	public void run(String[] args) throws Exception
	{
		//TODO: first pass, print list of variants that miss a CADD score
		//we can do SNVs easy but not indels and other 'non predictable' changes
		//send off to CADD webservice, and give the output so we use this in second pass
		
		if (!(args.length == 5))
		{
			throw new Exception("Must supply 5 arguments: mode ('"+MODE_ANALYSIS+"' or '"+MODE_CREATECADDFILE+"'), input VCF, patient groups, CCGG, CADD_suppl");
		}
		
		String mode = args[0];
		if(!mode.equals(MODE_ANALYSIS) && !mode.equals(MODE_CREATECADDFILE))
		{
			throw new Exception("Mode must be '"+MODE_ANALYSIS+"' or '"+MODE_CREATECADDFILE+"'");
		}

		File vcfFile = new File(args[1]);
		if (!vcfFile.isFile())
		{
			throw new Exception("Input VCF file does not exist or directory: " + vcfFile.getAbsolutePath());
		}

		File patientGroups = new File(args[2]);
		if (!patientGroups.isFile())
		{
			System.out.println("WARNING: patient groups file does not exist or directory: " + patientGroups.getAbsolutePath() +", not going to use patient groups!!");
		}
		
		File ccggFile = new File(args[3]);
		if (!ccggFile.isFile())
		{
			throw new Exception("CCGG file does not exist or directory: " + ccggFile.getAbsolutePath());
		}

		File caddFile = new File(args[4]);
		if(mode.equals(MODE_ANALYSIS))
		{
			if (!caddFile.isFile())
			{
				throw new Exception("CADD file does not exist or directory: " + ccggFile.getAbsolutePath());
			}
		}
		else if(mode.equals(MODE_CREATECADDFILE))
		{
			if (caddFile.isFile())
			{
				System.out.println("WARNING: going to create new CADD intermediate but already exists at " + caddFile.getAbsolutePath());
			}
		}
		
		FEEWAS cf = new FEEWAS();
		cf.go(mode, vcfFile, patientGroups, ccggFile, caddFile);

	}

}
