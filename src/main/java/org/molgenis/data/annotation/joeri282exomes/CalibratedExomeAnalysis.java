package org.molgenis.data.annotation.joeri282exomes;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Scanner;
import java.util.Set;

import org.molgenis.calibratecadd.support.LoadCADDWebserviceOutput;
import org.molgenis.data.Entity;
import org.molgenis.data.annotation.cmd.CommandLineAnnotatorConfig;
import org.molgenis.data.annotation.entity.impl.SnpEffAnnotator.Impact;
import org.molgenis.data.annotation.joeri282exomes.Judgment.Classification;
import org.molgenis.data.annotation.joeri282exomes.Judgment.Method;
import org.molgenis.data.vcf.VcfRepository;
import org.molgenis.data.vcf.utils.VcfUtils;
import org.springframework.context.annotation.AnnotationConfigApplicationContext;
import org.springframework.stereotype.Component;

@Component
public class CalibratedExomeAnalysis
{

	private static final String MODE_ANALYSIS = "analysis";
	private static final String MODE_CREATECADDFILE = "create_file_for_cadd";
	
	private Set<String> genesForVcfStats = new HashSet<String>();
	private Set<String> samplesForVcfStats = new HashSet<String>();
	
	private HashMap<String, String> sampleToGroup;
//	private HashMap<String, List<CandidateVariant>> sampleToPathogenicVariants;
	
	private List<CandidateVariant> pathoVariants;
	private List<CandidateVariant> vousVariants;
	private List<CandidateVariant> benignVariants;
	
	private int variantRefAltGenePatho = 0;
	private int variantRefAltGeneVOUS = 0;
	private int variantRefAltGeneBenign = 0;
	private int variantSkippedBadFilter = 0;
	private int noImpact = 0;
	private int noSamples = 0;
	private int totalPositionsSeen = 0;
	private int totalVariantRefAltGeneCombinationsSeen = 0;
	
	public void printResultsAsMatrix(List<CandidateVariant> pathoVariants, List<CandidateVariant> vousVariants, List<CandidateVariant> benignVariants)
	{
		System.out.print("Gene:variant");
		for(String sample : samplesForVcfStats)
		{
			System.out.print("\t" + sample);
		}
		System.out.print("\n");
		printVariantListAsMatrix(pathoVariants);
		printVariantListAsMatrix(vousVariants);
		printVariantListAsMatrix(benignVariants);
	}
	
	public void printResultsAsTable(List<CandidateVariant> pathoVariants, List<CandidateVariant> vousVariants, List<CandidateVariant> benignVariants) throws Exception
	{
		System.out.println("Chrom" + "\t" + "Pos" + "\t" + "Ref" + "\t" + "Alt" + "\t" + "AffectingAlt" + "\t" + "AltIndexInVCF" + "\t" + "AffectedGene" + "\t" + "cDNA" + "\t" + "aaChange" + "\t" + "Effect" + "\t" + "ID" + "\t" + "GoNL-AF" + "\t" + "ExAC-AF" + "\t" + "1000G-AF" + "\t" + "CADD-score" + "\t" +  "Classification" + "\t" + "Method" + "\t" + "Reason");
		printVariantListAsTable(pathoVariants);
		printVariantListAsTable(vousVariants);
		printVariantListAsTable(benignVariants);
	}
	
	public void printVariantListAsTable(List<CandidateVariant> variants) throws Exception
	{
		for(CandidateVariant cv : variants)
		{
			String chr = cv.getVcfRecord().getString("#CHROM");
			String pos = cv.getVcfRecord().getString("POS");
			String ref = cv.getVcfRecord().getString("REF");
			String alt = cv.getVcfRecord().getString("ALT");
			String altAff = cv.allele;
			int altIndex = cv.altIndex;
			String gene = cv.gene;
			Classification clsf = cv.judgment.classification;
			Method method = cv.judgment.method;
			String reason = cv.judgment.reason;
			String ann = CCGGUtils.getAnn(cv.getVcfRecord().getString("ANN"), cv.gene, cv.allele);
			String[] annSplit = ann.split("\\|", -1);
			String cDNA = annSplit[9];
			String aaChange = annSplit[10].isEmpty() ? "" : annSplit[10];
			String effect = annSplit[1];
			//String genomic = cv.vcfRecord.getString("#CHROM") + ":" + cv.vcfRecord.getString("POS") + " " + cv.vcfRecord.getString("REF") + " to " + cv.vcfRecord.getString("ALT");
			String id = (cv.vcfRecord.getString("ID") != null && !cv.vcfRecord.getString("ID").equals("") && !cv.vcfRecord.getString("ID").equals(".")) ? cv.vcfRecord.getString("ID") : "";
			String gonlAF = cv.vcfRecord.getString("GoNL_AF") == null ? "0" : (cv.vcfRecord.getString("GoNL_AF").split("\\|", -1)[altIndex].isEmpty() ? "0" : cv.vcfRecord.getString("GoNL_AF").split("\\|", -1)[altIndex]);
			String exacAF = cv.vcfRecord.getString("EXAC_AF") == null ? "0" : cv.vcfRecord.getString("EXAC_AF").replace(".,", "0,").replace(",.", ",0").split(",", -1)[altIndex];
			String _1kgAF = cv.vcfRecord.getString("Thousand_Genomes_AF") == null ? "0" : cv.vcfRecord.getString("Thousand_Genomes_AF").replace(".,", "0,").replace(",.", ",0").split(",", -1)[altIndex];
			Double cadd_STR = cv.cadd;
			System.out.println(chr + "\t" + pos + "\t" + ref + "\t" + alt + "\t" + altAff + "\t" + altIndex + "\t" + gene + "\t" + cDNA + "\t" + aaChange + "\t" + effect + "\t" + id + "\t" + gonlAF + "\t" + exacAF + "\t" + _1kgAF + "\t" + cadd_STR + "\t" + clsf + "\t" + method + "\t" + reason);
		}
	}
	
	public void printVariantListAsMatrix(List<CandidateVariant> variants)
	{
		for(CandidateVariant cv : variants)
		{
			String[] annSplit = cv.getVcfRecord().get("ANN").toString().split("\\|", -1);
			String cDNA = annSplit[9];
			System.out.print(cv.gene + ":" + cDNA);
			for(String sample : samplesForVcfStats)
			{
				if(cv.getSampleIds().containsKey(sample))
				{
					System.out.print("\t" + "1");
				}
				else
				{
					System.out.print("\t" + "0");
				}
				
			}
			System.out.print("\n");
		}
	}
	
	public void go(String mode, File vcfFile, File patientGroups, File ccggFile, File caddFile) throws Exception
	{
		long startTime = System.currentTimeMillis();
		pathoVariants = new ArrayList<CandidateVariant>();
		vousVariants = new ArrayList<CandidateVariant>();
		benignVariants = new ArrayList<CandidateVariant>();
		loadSampleToGroup(patientGroups);

		CCGGUtils ccgg = new CCGGUtils(ccggFile);
		
		BufferedWriter pathogenicVariantsVCF = new BufferedWriter(new PrintWriter(new File(vcfFile.getParent(), "pathogenicVariants.vcf")));
		BufferedWriter benignVariantsVCF = new BufferedWriter(new PrintWriter(new File(vcfFile.getParent(), "benignVariants.vcf")));
		BufferedWriter vousVariantsVCF = new BufferedWriter(new PrintWriter(new File(vcfFile.getParent(), "vousVariants.vcf")));
		
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
			if(totalPositionsSeen % 1000 == 0)
			{
				System.out.println("Seen " + totalPositionsSeen + " variants, in a total of " + totalVariantRefAltGeneCombinationsSeen + " ref/alt/gene combinations. Classified " + variantRefAltGenePatho+ " as pathogenic, " + variantRefAltGeneBenign + " as benign, failed to classify/VOUS " + (noSamples+noImpact+variantRefAltGeneVOUS) + ".");
			}
			totalPositionsSeen++;
			
			Entity record = it.next();

			String filter = record.getString("FILTER");

			if (filter != null && !filter.equals("PASS"))
			{
				variantSkippedBadFilter++;
				continue;
			}
			
		//	System.out.println(record.toString());
			String chr = record.getString("#CHROM");
			String pos = record.getString("POS");
			String ref = record.getString("REF");
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
			
			Set<String> genes = CCGGUtils.getGenesFromAnn(ann);
			
			
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
				
				ArrayList<CandidateVariant> multipleCandidatesForAllele = new ArrayList<CandidateVariant>();
				for(String gene : genes)
				{
					genesForVcfStats.add(gene);
					totalVariantRefAltGeneCombinationsSeen++;
					
					Impact impact = CCGGUtils.getImpact(ann, gene, alt);
					if(impact == null)
					{
						noImpact++;
						continue;
					}
					
					Judgment judgment = null;
				
					judgment = ccgg.classifyVariant(gene, exac_af, impact, cadd);

					if (judgment.classification.equals(Judgment.Classification.Pathogn) || judgment.classification.equals(Judgment.Classification.VOUS))
					{
						HashMap<String, Entity> samples = findInterestingSamples(record, i, exac_het, exac_hom);
						if(samples.size() > 0)
						{
							CandidateVariant cv = new CandidateVariant(record, alt, i, gene, cadd, judgment, samples);
							if(judgment.classification.equals(Judgment.Classification.Pathogn))
							{
								pathoVariants.add(cv);
								variantRefAltGenePatho++;
								VcfUtils.writeToVcf(record, pathogenicVariantsVCF);
								pathogenicVariantsVCF.write('\n');
							}
							if(judgment.classification.equals(Judgment.Classification.VOUS))
							{
								vousVariants.add(cv);
								variantRefAltGeneVOUS++;
								VcfUtils.writeToVcf(record, vousVariantsVCF);
								vousVariantsVCF.write('\n');
							}
						}
						else
						{
							noSamples++;
						}
					}
					else if(judgment.classification.equals(Judgment.Classification.Benign))
					{
						HashMap<String, Entity> samples = findInterestingSamples(record, i, exac_het, exac_hom);
						CandidateVariant cv = new CandidateVariant(record, alt, i, gene, cadd, judgment, samples);
						benignVariants.add(cv);
						variantRefAltGeneBenign++;
						VcfUtils.writeToVcf(record, benignVariantsVCF);
						benignVariantsVCF.write('\n');
					}
				}
				
			}
		}
		vcf.close();
		pathogenicVariantsVCF.close();
		benignVariantsVCF.close();
		vousVariantsVCF.close();
		
		if(mode.equals(MODE_CREATECADDFILE))
		{
			pw.flush();
			pw.close();
		}

	//	printResults(pathoVariants, "Pathogenic");
	//	printResults(vousVariants, "VOUS");
		
	//	printResultsAsMatrix(pathoVariants, vousVariants, benignVariants);
		
		printResultsAsTable(pathoVariants, vousVariants, benignVariants);
		
		long endTime = System.currentTimeMillis();
		
		System.out.println("\n## COUNTS ##\n");
		System.out.println("total nr of samples in vcf = " + samplesForVcfStats.size());
		System.out.println("total nr of genes in vcf = " + genesForVcfStats.size());
		System.out.println("total positions seen = " + totalPositionsSeen);
		System.out.println("total positions that did not pass QC and were skipped: " + variantSkippedBadFilter);
		System.out.println("total variant combinations of alleles and genes seen = " + totalVariantRefAltGeneCombinationsSeen);
		System.out.println("any combinations classified benign = " + variantRefAltGeneBenign);
		System.out.println("any combinations classified pathogenic = " + variantRefAltGenePatho);
		System.out.println("classification failed / VOUS = " + variantRefAltGeneVOUS);
		System.out.println("classification failed (no samples) = " + noSamples);
		System.out.println("classification failed (no impact) = " + noImpact);
		System.out.println("runtime: " + (endTime-startTime)/1000.0 + " seconds, " + totalVariantRefAltGeneCombinationsSeen/((endTime-startTime)/1000.0) + " interpretations per second");
	
	}
	
	public void printResults(List<CandidateVariant> candVariants, String name) throws Exception
	{
		System.out.println("\n## RESULTS ##\n");
		
		for(Method conf : Method.values())
		{
			System.out.println("# " + name + " candidates, method: " + conf + "\n");
			for(CandidateVariant cv : candVariants)
			{
				if(cv.judgment.method.equals(conf))
				{
			//		System.out.println(cv.toString());
					
					String ann = CCGGUtils.getAnn(cv.getVcfRecord().getString("ANN"), cv.gene, cv.allele);
					String[] annSplit = ann.split("\\|", -1);
					String cDNA = annSplit[9];
					String aaChange = annSplit[10].isEmpty() ? null : annSplit[10];
					String effect = annSplit[1];
					String genomic = cv.vcfRecord.getString("#CHROM") + ":" + cv.vcfRecord.getString("POS") + " " + cv.vcfRecord.getString("REF") + " to " + cv.vcfRecord.getString("ALT");
					String id = (cv.vcfRecord.getString("ID") != null && !cv.vcfRecord.getString("ID").equals("") && !cv.vcfRecord.getString("ID").equals(".")) ? cv.vcfRecord.getString("ID") : null;
					String gonlAF = cv.vcfRecord.getString("GoNL_AF") == null ? "0" : cv.vcfRecord.getString("GoNL_AF");
					String exacAF = cv.vcfRecord.getString("EXAC_AF") == null ? "0" : cv.vcfRecord.getString("EXAC_AF");
					String _1kgAF = cv.vcfRecord.getString("Thousand_Genomes_AF") == null ? "0" : cv.vcfRecord.getString("Thousand_Genomes_AF");
					
					System.out.println("Variant: " + cv.gene + ":" + cDNA + (aaChange != null || id != null ? " (" : "") + (aaChange != null ? aaChange + (id != null ? ", " + id: "") : "") + (aaChange != null || id != null ? ")" : "") + ", genomic: " +genomic + " (pathogenic allele: " + cv.allele + ")");
					System.out.println("Effect: " + effect + ", GoNL MAF: " + gonlAF+ ", ExAC MAF: " + exacAF+ ", 1KG MAF: " + _1kgAF);
					System.out.println("Reason: " + cv.getJudgment().reason);
					
					System.out.println("Samples carrying this allele:");
					for(String sampleId : cv.getSampleIds().keySet())
					{
						System.out.println(sampleId + ", genotype: " + cv.getSampleIds().get(sampleId).getString("GT") + ", cov: " + cv.getSampleIds().get(sampleId).getString("DP") + ", allelic cov: " + cv.getSampleIds().get(sampleId).getString("AD"));
					}
					System.out.println();
				}
				
			}
		}
	}
	
	public HashMap<String, Entity> findInterestingSamples(Entity record, int altIndex, int exacHets, int exacHoms)
	{
		HashMap<String, Entity> homSampleIds = new HashMap<String, Entity>();
		HashMap<String, Entity> hetSampleIds = new HashMap<String, Entity>();
		Iterable<Entity> sampleEntities = record.getEntities(VcfRepository.SAMPLES);
		//first alt has index 0, but we want to check 0/1, 1/1 etc. So +1.
		altIndex = altIndex + 1;
		
		int hetCount = 0;
		int homCount = 0;
		
		for (Entity sample : sampleEntities)
		{
			String genotype = sample.get("GT").toString();
			
			String sampleName = sample.get("ORIGINAL_NAME").toString();
			samplesForVcfStats.add(sampleName);

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

			
			if (genotype.equals("0/" + altIndex) || genotype.equals(altIndex + "/0")
					|| genotype.equals("0|" + altIndex) || genotype.equals(altIndex + "|0") )
			{
				//interesting!
				
				hetSampleIds.put(sampleName, sample);
				hetCount++;
			}
			
			
			if ( genotype.equals(altIndex + "/" + altIndex) || genotype.equals(altIndex + "|" + altIndex) )
			{
				//interesting!
				homSampleIds.put(sampleName, sample);
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
		HashMap<String, Entity> sampleIds = new HashMap<String, Entity>();
		sampleIds.putAll(homSampleIds);
		sampleIds.putAll(hetSampleIds);
		return sampleIds;
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
			throw new Exception("Patient groups file does not exist or directory: " + patientGroups.getAbsolutePath());
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
		
		CalibratedExomeAnalysis cf = new CalibratedExomeAnalysis();
		cf.go(mode, vcfFile, patientGroups, ccggFile, caddFile);

	}

}
