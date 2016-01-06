package org.molgenis.data.annotation.joeri282exomes;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Scanner;
import java.util.Set;

import org.molgenis.calibratecadd.support.LoadCADDWebserviceOutput;
import org.molgenis.calibratecadd.support.VariantClassificationException;
import org.molgenis.data.Entity;
import org.molgenis.data.annotation.cmd.CommandLineAnnotatorConfig;
import org.molgenis.data.annotation.entity.impl.SnpEffAnnotator.Impact;
import org.molgenis.data.annotation.joeri282exomes.Judgment.Classification;
import org.molgenis.data.annotation.joeri282exomes.Judgment.Method;
import org.molgenis.data.vcf.VcfRepository;
import org.springframework.context.annotation.AnnotationConfigApplicationContext;
import org.springframework.stereotype.Component;

@Component
public class CalibratedExomeAnalysis
{

	private static final String MODE_ANALYSIS = "analysis";
	private static final String MODE_CREATECADDFILE = "create_file_for_cadd";
	
	private HashMap<String, String> sampleToGroup;
//	private HashMap<String, List<CandidateVariant>> sampleToPathogenicVariants;
	
	private List<CandidateVariant> pathogenicVariants;
	private List<Entity> noJudgmentVariants;
	private List<Entity> conflictingJudgmentVariants;
	
	private int variantRefAltGenePatho = 0;
	private int variantRefAltGeneBenign = 0;
	private int variantSkippedBadFilter = 0;
	private int failedToClassify = 0;
	private int conflictingClassifications = 0;
	private int totalVariantSeen = 0;
	private int totalVariantRefAltGeneCombinationsSeen = 0;
	
	public void go(String mode, File vcfFile, File patientGroups, File ccggFile, File caddFile) throws Exception
	{
		pathogenicVariants = new ArrayList<CandidateVariant>();
		noJudgmentVariants = new ArrayList<Entity>();
		conflictingJudgmentVariants = new ArrayList<Entity>();
		
		loadSampleToGroup(patientGroups);

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
				System.out.println("Seen " + totalVariantSeen + " variants, in a total of " + totalVariantRefAltGeneCombinationsSeen + " ref/alt/gene combinations. Classified " + variantRefAltGenePatho+ " as pathogenic, " + variantRefAltGeneBenign + " as benign, failed to classify " + failedToClassify + ".");
			}
			totalVariantSeen++;
			
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
						pw.println(chr + "\t" + pos + "\t" + "." + "\t" + ref + "\t" + alt);
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
							System.out.println("WARNING: CADD score missing for " + chr + " " + pos + " " + ref + " " + alt + " !");
						}
					}
				}
				
				ArrayList<CandidateVariant> multipleCandidatesForAllele = new ArrayList<CandidateVariant>();
				for(String gene : genes)
				{
					totalVariantRefAltGeneCombinationsSeen++;
					
					Impact impact = CCGGUtils.getImpact(ann, gene, alt);
					Judgment judgment = null;
					try
					{
						judgment = ccgg.classifyVariant(gene, exac_af, impact, cadd);

						HashMap<String, Entity> samples = findInterestingSamples(record, i, exac_het, exac_hom);
						if (samples.size() > 0)
						{
							CandidateVariant cv = new CandidateVariant(record, alt, gene, judgment, samples);
							multipleCandidatesForAllele.add(cv);
						}

					}
					catch(VariantClassificationException e)
					{
						
					}

				}
				
				//if no judgment, add null for this variant
				if(multipleCandidatesForAllele.size() == 0)
				{
					System.out.println("WARNING: no classification could be made for " + chr + ":" + pos + " " + ref + "/" + alt);
					failedToClassify++;
					noJudgmentVariants.add(record);
					continue;
				}
				
				
				//go through the possible classifications and check if any of them are conflicting
				//also, if we have a calibrated judgment, 
				int nrOfBenignClsf = 0;
				int nrOfPathognClsf = 0;
				boolean hasCalibratedJudgment = false;
				for(CandidateVariant cv : multipleCandidatesForAllele)
				{
					if(cv.judgment.getClassification().equals(Classification.Benign))
					{
						nrOfBenignClsf++;
					}
					if(cv.judgment.getClassification().equals(Classification.Pathogn))
					{
						nrOfPathognClsf++;
					}
					if(cv.judgment.getConfidence().equals(Method.calibrated))
					{
						hasCalibratedJudgment = true;
					}
				}
				
				//check if we have any conflicts
				if(nrOfBenignClsf > 0 && nrOfPathognClsf > 0)
				{
					System.out.println("WARNING: conflicting classification! adding no judgment for this variant: " + chr + ":" + pos + " " + ref + "/" + alt + ", judged: " + multipleCandidatesForAllele.toString() );
					conflictingClassifications++;
					conflictingJudgmentVariants.add(record);
				}
				else
				{
					for(CandidateVariant cv : multipleCandidatesForAllele)
					{
						if(cv.judgment.classification.equals(Classification.Benign))
						{
							variantRefAltGeneBenign++;
							break;
						}
						//if we know we have calibrated results, wait for it, then add it, and then break
						if(hasCalibratedJudgment && cv.judgment.getConfidence().equals(Method.calibrated) && cv.judgment.classification.equals(Classification.Pathogn))
						{
							variantRefAltGenePatho++;
							pathogenicVariants.add(cv);
							break;
						}
						//if not, might as well add this one and be done
						else if(!hasCalibratedJudgment && cv.judgment.classification.equals(Classification.Pathogn))
						{
							variantRefAltGenePatho++;
							pathogenicVariants.add(cv);
							break;
						}
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

		printResults(pathogenicVariants);
		
		System.out.println("\n## COUNTS ##\n");
		System.out.println("total variant seen = " + totalVariantSeen);
		System.out.println("variants that did not pass QC and were skipped: " + variantSkippedBadFilter);
		System.out.println("total variants incl. altallele/gene combinations seen = " + totalVariantRefAltGeneCombinationsSeen);
		System.out.println("relevant combinations, classified benign = " + variantRefAltGeneBenign);
		System.out.println("relevant combinations, classified pathogenic = " + variantRefAltGenePatho);
		System.out.println("classification failed = " + failedToClassify);
		System.out.println("classification conflicting = " + conflictingClassifications);

	}
	
	public void printResults(List<CandidateVariant> pathogenicVariants) throws Exception
	{
		System.out.println("\n## RESULTS ##\n");
		
		for(Method conf : Method.values())
		{
			System.out.println("\nPathogenic candidates, method: " + conf);
			for(CandidateVariant cv : pathogenicVariants)
			{
				if(cv.judgment.method.equals(conf))
				{
			//		System.out.println(cv.toString());
					
					String ann = CCGGUtils.getAnn(cv.getVcfRecord().getString("ANN"), cv.gene, cv.allele);
					String[] annSplit = ann.split("\\|", -1);
					String cDNA = annSplit[9];
					String aaChange = annSplit[10];
					String effect = annSplit[1];
					String genomic = cv.vcfRecord.getString("#CHROM") + ":" + cv.vcfRecord.getString("POS") + " " + cv.vcfRecord.getString("REF") + " to " + cv.vcfRecord.getString("ALT");
					String id = (cv.vcfRecord.getString("ID") != null && !cv.vcfRecord.getString("ID").equals("") && !cv.vcfRecord.getString("ID").equals(".")) ? (", " + cv.vcfRecord.getString("ID")) : "";
					String gonlAF = cv.vcfRecord.getString("GoNL_AF").equals(".") ? "0" : cv.vcfRecord.getString("GoNL_AF");
					String exacAF = cv.vcfRecord.getString("EXAC_AF") == null ? "0" : cv.vcfRecord.getString("EXAC_AF");
					String _1kgAF = cv.vcfRecord.getString("Thousand_Genomes_AF") == null ? "0" : cv.vcfRecord.getString("Thousand_Genomes_AF");
					
					System.out.println("Variant: " + cv.gene + ":" + cDNA + " (" + aaChange + id + "), genomic: " +genomic + " (pathogenic allele: " + cv.allele + ")");
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
				hetSampleIds.put(sampleName, sample);
				hetCount++;
			}
			
			
			if ( genotype.equals(altIndex + "/" + altIndex) || genotype.equals(altIndex + "|" + altIndex) )
			{
				//interesting!
				String sampleName = sample.get("ORIGINAL_NAME").toString();
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
