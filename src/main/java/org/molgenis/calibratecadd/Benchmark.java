package org.molgenis.calibratecadd;

import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Set;

import org.apache.commons.lang3.EnumUtils;
import org.molgenis.calibratecadd.support.*;
import org.molgenis.calibratecadd.support.JudgedVariant.ExpertClassification;
import org.molgenis.data.Entity;

import org.molgenis.data.annotation.entity.impl.gavin.GavinAlgorithm;
import org.molgenis.data.annotation.entity.impl.gavin.GavinAnnotator;
import org.molgenis.data.annotation.entity.impl.gavin.GavinEntry;
import org.molgenis.data.annotation.entity.impl.gavin.Judgment;
import org.molgenis.data.annotation.entity.impl.snpEff.Impact;
import org.molgenis.data.vcf.VcfRepository;

import org.molgenis.data.annotation.entity.impl.gavin.Judgment.Classification;
import org.molgenis.data.annotation.entity.impl.gavin.Judgment.Method;

/**
* Assess the performance of an in silico variant pathogenicity prediction tool on gold standard datasets.
* 
* In short: you can run any 'gold standard file' against the CCGG tool, as long as this is a VCF with columns for MVL, CLSF, CADD_SCALED and EXAC_AF.
* If you want to check the performance of a different tool, you must score the MVL manually and link to the (minimized/cleaned) output of these tools here.
* This emulates the classification process using a different tool.
* You can do this already for the GoldStandard_UMCG_MVLs_noClinVar data set. Use the "Prediction_UMCG_MVLs_noClinVar_XXXX" files for this (see code below).
* 
* First argument: the GAVIN output file. Typically: GAVIN_calibrations_r0.2.tsv
* 
* Second argument:
* A gold standard file to test against. Typically a VCF file with 'MVL' and 'CLSF' columns denoting expert validated interpretations.
* Three are used and available here:
* 1) MutationTaster benchmark set (more info: http://www.mutationtaster.org/info/Comparison_20130328_with_results_ClinVar.html)
* 		GoldStandard_MutationTaster_TestSet.vcf
* 2) VariBench benign and pathogenic datasets. Used to train and test PON-P2.
* 		GoldStandard_VariBench_ToleranceDS7_all.vcf.gz
* 3) UMCG Managed Variant Lists for 5 disease types.
* 		GoldStandard_UMCG_MVLs_noClinVar.vcf
* 
* Third argument:
* The tool to classify the variants with and subsequently report the performance of this method.
* Available by default are: 'ccgg' and 'naive', using CCGG_ClassificationSource_GeneSummary.tsv
* Requiring extra input (pre-scored variants in minimized files) are: 'ponp2', 'mutationtaster2', 'provean', 'sift'
* See file examples ('Prediction_XXX') in the hardcoded locations on how to plug in different files for different MVL inputs.
* 
* Example usages:
* 
 /Users/jvelde/github/maven/gavin/data/CCGG_ClassificationSource_GeneSummary.tsv
 /Users/jvelde/github/maven/gavin/data/GoldStandard_UMCG_MVLs_noClinVar.vcf
 ccgg
*
 /Users/jvelde/github/maven/gavin/data/CCGG_ClassificationSource_GeneSummary.tsv
 /Users/jvelde/github/maven/gavin/data/GoldStandard_MutationTaster_TestSet.vcf
 ccgg
*
 /Users/jvelde/github/maven/gavin/data/CCGG_ClassificationSource_GeneSummary.tsv
 /Users/jvelde/github/maven/gavin/data/GoldStandard_UMCG_MVLs_noClinVar.vcf
 provean
*
 /Users/jvelde/github/maven/gavin/data/CCGG_ClassificationSource_GeneSummary.tsv
 /Users/jvelde/github/maven/gavin/data/GoldStandard_VariBench_ToleranceDS7_all.vcf.gz
 ccgg
* 
 /Users/jvelde/github/maven/gavin/data/CCGG_ClassificationSource_GeneSummary.tsv
 /Users/jvelde/github/maven/gavin/data/GoldStandard_UMCG_MVLs_noClinVar.vcf
 ponp2
* 
* MS Windows paths:
 C:\Users\Joeri\github\gavin\data\CCGG_ClassificationSource_GeneSummary.tsv
 C:\Users\Joeri\github\gavin\data\GoldStandard_UMCG_MVLs_noClinVar.vcf
 ccgg
* 
* @author jvelde
*
*/
public class Benchmark
{
	public enum ToolNames{
		GAVIN, GAVINnocal, PONP2, CADD_Thr15, CADD_Thr20, CADD_Thr25, PROVEAN, SIFT, PolyPhen2, MSC_ClinVar95CI, MSC_HGMD99CI, Condel, PredictSNP2, FATHMM, GWAVA, FunSeq, DANN
	}
	
	public static void main(String[] args) throws Exception
	{
		if(args.length != 4)
		{
			throw new Exception("please provide: gavin file, variant vcf, tool, output file, version [e.g. 'r0.2']");
		}
		new File(args[3]).createNewFile();
		new Benchmark(args[0], args[1], args[2], args[3], args[4]);
		GavinAnnotator ga = new GavinAnnotator();
	//	ga.gavin().annotateEntity();
	}

	HashMap<String, GavinEntry> gavinData;
	GavinAlgorithm gavin;
	HashMap<String, List<JudgedVariant>> judgedMVLVariants = new HashMap<String, List<JudgedVariant>>();
	int judgmentsInCalibratedGenes = 0;
	
	/**
	 * Take MVL (annotated with CADD, ExAC, SnpEff)
	 * Take CCGG thresholds
	 * Check if classification matches
	 * @throws Exception 
	 */
	public Benchmark(String predictionToolPath, String mvlLoc, String mode, String outFile, String version) throws Exception
	{
		if(!EnumUtils.isValidEnum(ToolNames.class, mode))
		{
			throw new Exception("mode needs to be one of : " + java.util.Arrays.asList(ToolNames.values()));
		}
		File mvlFile = new File(mvlLoc);
		if(!mvlFile.exists())
		{
			throw new Exception("MVL file "+mvlFile+" does not exist or is directory");
		}
		this.gavinData = loadCCGG(predictionToolPath + File.separatorChar + "GAVIN_calibrations_"+version+".tsv").getGeneToEntry();
		this.gavin = new GavinAlgorithm();
		scanMVL(mvlFile, predictionToolPath, ToolNames.valueOf(mode));
		ProcessJudgedVariantMVLResults.printResults(judgedMVLVariants, mode, mvlFile.getName(), judgmentsInCalibratedGenes, outFile);
	}
	
	public void scanMVL(File mvlFile, String predictionToolPath, ToolNames mode) throws Exception
	{
		
		VcfRepository vcfRepo = new VcfRepository(mvlFile, "mvl");
		
		java.util.Iterator<Entity> vcfRepoIter = vcfRepo.iterator();
		
		/**
		 * Here, we let PON-P2 or other tools handle the classification and see what happens
		 * Obviously, before this works, you must score your variant list with PON-P2 or other tools
		 * And then link this file below
		 */
		PONP2Results p2r = null;
		MutationTaster2Results m2r = null;
		ProveanAndSiftResults ps2r = null;
		PolyPhen2Results pf2r = null;
		MSCResults mscr = null;
		CondelResults condelr = null;
		PredictSNP2Results predictSNP2r = null;
		if (mode.equals(ToolNames.PONP2))
		{
			p2r = new PONP2Results(new File(predictionToolPath, "PON-P2.tsv"));
		}
		if (mode.equals(ToolNames.PROVEAN) || mode.equals(ToolNames.SIFT))
		{
			ps2r = new ProveanAndSiftResults(new File(predictionToolPath, "PROVEAN_SIFT.tsv"));
		}
		if (mode.equals(ToolNames.PolyPhen2))
		{
			pf2r = new PolyPhen2Results(new File(predictionToolPath, "PolyPhen2.tsv"));
		}
		if (mode.equals(ToolNames.MSC_ClinVar95CI))
		{
			mscr = new MSCResults(new File(predictionToolPath, "MSC_CADD_cutoffs_ClinVar95CI.tsv"));
		}
		if (mode.equals(ToolNames.MSC_HGMD99CI))
		{
			mscr = new MSCResults(new File(predictionToolPath, "MSC_CADD_cutoffs_HGMD99CI.tsv"));
		}
		if (mode.equals(ToolNames.Condel))
		{
			condelr = new CondelResults(new File(predictionToolPath, "Condel.tsv"));
		}
		if (mode.equals(ToolNames.PredictSNP2) || mode.equals(ToolNames.FATHMM) || mode.equals(ToolNames.FunSeq) || mode.equals(ToolNames.GWAVA) || mode.equals(ToolNames.DANN) )
		{
			predictSNP2r = new PredictSNP2Results(new File(predictionToolPath, "PredictSNP2.vcf.gz"));
		}
		
		System.out.println("Running in mode: " + mode);

		while (vcfRepoIter.hasNext())
		{
			Entity record = vcfRepoIter.next();
			
			String chr = record.getString("#CHROM");
			String pos = record.getString("POS");
			String ref = record.getString("REF");
			String alt = record.getString("ALT");
			if(alt.contains(","))
			{
				throw new Exception("Did not expect multiple alt alleles! " + record.toString());
			}
			
			Double getMAF = GavinUtils.getInfoForAllele(record, "EXAC_AF", alt);
			double MAF = getMAF == null ? 0 : getMAF;
			Double CADDscore = GavinUtils.getInfoForAllele(record, "CADD_SCALED", alt);
			String ann = record.getString("ANN");
			Set<String> genes = GavinUtils.getGenesFromAnn(ann);
			String id = record.getString("ID");
			
			//for some variants, we have GENE:MUTATION in the ID
			// for example, in the MVL data, we see "GUSB:c.1222C>T", and in the VariBench data: "PAH:c.617A>G"
			//if this is present, we use this to narrow the scope by matching the annotated genes to the gene symbol here
			String[] idSplit = id.split(":", -1);
			boolean hasGeneId = idSplit.length > 1 ? true : false;
			String geneFromId = idSplit[0];
			
			String mvlClassfc = record.getString("CLSF");
			String mvlName = record.getString("MVL");
			
			ArrayList<Judgment> multipleJudgments = new ArrayList<Judgment>();
			
			boolean geneToIdMatchFound = false;
			for(String gene : genes)
			{
				if(!hasGeneId || gene.equals(geneFromId))
				{
					geneToIdMatchFound = true;
					Impact impact = GavinUtils.getImpact(ann, gene, alt);

					Judgment judgment;
					if (mode.equals(ToolNames.GAVIN))
					{
						judgment = gavin.classifyVariant(impact, CADDscore, MAF, gene, null, this.gavinData);
					}
					else if(mode.equals(ToolNames.GAVINnocal))
					{
						judgment = gavin.genomewideClassifyVariant(impact, CADDscore, MAF, gene);
					}
					else if (mode.equals(ToolNames.PONP2))
					{
						judgment = p2r.classifyVariantUsingPONP2Results(chr, pos, ref, alt);
					}
					else if (mode.equals(ToolNames.PROVEAN))
					{
						judgment = ps2r.classifyVariantUsingProveanResults(chr, pos, ref, alt);
					}
					else if (mode.equals(ToolNames.SIFT))
					{
						judgment = ps2r.classifyVariantUsingSiftResults(chr, pos, ref, alt);
					}
					else if (mode.equals(ToolNames.PolyPhen2))
					{
						judgment = pf2r.classifyVariantUsingPolyPhen2Results(chr, pos, ref, alt);
					}
					else if (mode.equals(ToolNames.MSC_ClinVar95CI) || mode.equals(ToolNames.MSC_HGMD99CI))
					{
						judgment = mscr.classifyVariantUsingMSCResults(gene, CADDscore);
					}
					else if (mode.equals(ToolNames.Condel))
					{
						judgment = condelr.classifyVariantUsingCondelResults(chr, pos, ref, alt);
					}
					else if (mode.equals(ToolNames.PredictSNP2) || mode.equals(ToolNames.FATHMM) || mode.equals(ToolNames.FunSeq) || mode.equals(ToolNames.GWAVA) || mode.equals(ToolNames.DANN) )
					{
						judgment = predictSNP2r.classifyVariantUsingPredictSNP2Results(chr, pos, ref, alt, mode.toString());
					}
					else if (mode.equals(ToolNames.CADD_Thr15))
					{
						if(CADDscore != null && CADDscore > 15)
						{
							judgment = new Judgment(Classification.Pathogenic, Method.calibrated, gene, "CADD score > 15");
						}
						else if(CADDscore != null && CADDscore <= 15)
						{
							judgment = new Judgment(Classification.Benign, Method.calibrated, gene, "CADD score <= 15");
						}
						else{
							judgment = new Judgment(Classification.VOUS, Method.calibrated, gene, "CADD score not available");
						}
					}
					else if (mode.equals(ToolNames.CADD_Thr20))
					{
						if(CADDscore != null && CADDscore > 20)
						{
							judgment = new Judgment(Classification.Pathogenic, Method.calibrated, gene, "CADD score > 20");
						}
						else if(CADDscore != null && CADDscore <= 20)
						{
							judgment = new Judgment(Classification.Benign, Method.calibrated, gene, "CADD score <= 20");
						}
						else{
							judgment = new Judgment(Classification.VOUS, Method.calibrated, gene, "CADD score not available");
						}
					}
					else if (mode.equals(ToolNames.CADD_Thr25))
					{
						if(CADDscore != null && CADDscore > 25)
						{
							judgment = new Judgment(Classification.Pathogenic, Method.calibrated, gene, "CADD score > 25");
						}
						else if(CADDscore != null && CADDscore <= 25)
						{
							judgment = new Judgment(Classification.Benign, Method.calibrated, gene, "CADD score <= 25");
						}
						else{
							judgment = new Judgment(Classification.VOUS, Method.calibrated, gene, "CADD score not available");
						}
					}
					else
					{
						throw new Exception("Mode unknown: " + mode);
					}
					
					multipleJudgments.add(judgment);
				}
			}
			if(hasGeneId && !geneToIdMatchFound)
			{
				if(mode.equals(ToolNames.GAVIN) && gavinData.containsKey(geneFromId)) { judgmentsInCalibratedGenes++; }
				System.out.println("WARNING: bad data for variant " + chr + ":" + pos + " " + ref + "/" + alt + ", no match from ID field gene to snpeff annotations!");
				multipleJudgments.add(new Judgment(Classification.VOUS, Method.calibrated, geneFromId, "Bad data!"));
			}
			
			//if no judgment, add null for this variant
			if(multipleJudgments.size() == 0)
			{
				throw new Exception("No judgments! should not occur.");
			}
			
			//go through the possible classifications and check if any of them are conflicting
			//also, if we have a calibrated judgment, 
			int nrOfBenignClsf = 0;
			int nrOfPathognClsf = 0;
			boolean hasCalibratedJudgment = false;
			for(Judgment judgment : multipleJudgments)
			{
				if(judgment.getClassification().equals(Classification.Benign))
				{
					nrOfBenignClsf++;
				}
				if(judgment.getClassification().equals(Classification.Pathogenic))
				{
					nrOfPathognClsf++;
				}
				if(judgment.getConfidence().equals(Method.calibrated))
				{
					hasCalibratedJudgment = true;
				}
			}
			
			/**
			 * Now we can assign the final verdict for this variant
			 */
			
			//check if we have any conflicts
			//TODO could be improved by prioritizing calibrated over genomewide results for our method
			if(nrOfBenignClsf > 0 && nrOfPathognClsf > 0)
			{
				if(mode.equals(ToolNames.GAVIN) && hasCalibratedJudgment) { judgmentsInCalibratedGenes++; }
				System.out.println("WARNING: conflicting classification! adding no judgment for this variant: " + chr + ":" + pos + " " + ref + "/" + alt + ", judged: " + multipleJudgments.toString() );
				addToMVLResults(new Judgment(Classification.VOUS, (hasCalibratedJudgment ? Method.calibrated : Method.genomewide), geneFromId, "Conflicting classification!!"), mvlClassfc, mvlName, record);
			}
			else
			{
				for(Judgment judgment : multipleJudgments)
				{
					//if we know we have calibrated results, wait for it, then add it, and then break
					if(hasCalibratedJudgment && judgment.getConfidence().equals(Method.calibrated))
					{
						addToMVLResults(judgment, mvlClassfc, mvlName, record);
						if(mode.equals(ToolNames.GAVIN)) { judgmentsInCalibratedGenes++; }
						break;
					}
					//if not, might as well add this one and be done
					//TODO: this means there may be multiple verdicts, e.g. 2x BENIGN for context in two genes, but we only add 1 of them, to keep things a bit more simple
					else if(!hasCalibratedJudgment)
					{
						addToMVLResults(judgment, mvlClassfc, mvlName, record);
						break;
					}
				}
			}
		}
	}
	
	private void addToMVLResults(Judgment judgment, String mvlClasfc, String mvl, Entity variant)
	{
		if(!judgedMVLVariants.containsKey(mvl))
		{
			List<JudgedVariant> variants = new ArrayList<JudgedVariant>();
			judgedMVLVariants.put(mvl, variants);
		}
		List<JudgedVariant> variants = judgedMVLVariants.get(mvl);
		variants.add(new JudgedVariant(judgment, variant, ExpertClassification.valueOf(mvlClasfc)));
	}
	
	public static GavinUtils loadCCGG(String ccggLoc) throws Exception
	{
		File ccggFile = new File(ccggLoc);
		if(!ccggFile.exists())
		{
			throw new Exception("CCGG file "+ccggFile+" does not exist or is directory");
		}
		GavinUtils ccggUtils = new GavinUtils(ccggFile);
		return ccggUtils;
	}
}
