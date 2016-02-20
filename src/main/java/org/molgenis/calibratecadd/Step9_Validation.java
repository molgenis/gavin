package org.molgenis.calibratecadd;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Set;

import org.molgenis.calibratecadd.support.CondelResults;
import org.molgenis.calibratecadd.support.MSCResults;
import org.molgenis.calibratecadd.support.MutationTaster2Results;
import org.molgenis.calibratecadd.support.PONP2Results;
import org.molgenis.calibratecadd.support.PolyPhen2Results;
import org.molgenis.calibratecadd.support.ProveanAndSiftResults;
import org.molgenis.calibratecadd.support.VariantClassificationException;
import org.molgenis.calibratecadd.support.JudgedVariant;
import org.molgenis.calibratecadd.support.JudgedVariant.ExpertClassification;
import org.molgenis.calibratecadd.support.ProcessJudgedVariantMVLResults;
import org.molgenis.data.Entity;
import org.molgenis.data.annotation.entity.impl.SnpEffAnnotator.Impact;
import org.molgenis.data.annotation.joeri282exomes.CCGGUtils;
import org.molgenis.data.annotation.joeri282exomes.Judgment;
import org.molgenis.data.annotation.joeri282exomes.Judgment.Classification;
import org.molgenis.data.annotation.joeri282exomes.Judgment.Method;
import org.molgenis.data.vcf.VcfRepository;

/**
* Assess the performance of an in silico variant pathogenicity prediction tool on gold standard datasets.
* 
* In short: you can run any 'gold standard file' against the CCGG tool, as long as this is a VCF with columns for MVL, CLSF, CADD_SCALED and EXAC_AF.
* If you want to check the performance of a different tool, you must score the MVL manually and link to the (minimized/cleaned) output of these tools here.
* This emulates the classification process using a different tool.
* You can do this already for the GoldStandard_UMCG_MVLs_noClinVar data set. Use the "Prediction_UMCG_MVLs_noClinVar_XXXX" files for this (see code below).
* 
* First argument: the CCGG output file. Typically: CCGG_ClassificationSource_GeneSummary.tsv
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
 /Users/jvelde/github/maven/molgenis-data-cadd/data/CCGG_ClassificationSource_GeneSummary.tsv
 /Users/jvelde/github/maven/molgenis-data-cadd/data/GoldStandard_UMCG_MVLs_noClinVar.vcf
 ccgg
*
 /Users/jvelde/github/maven/molgenis-data-cadd/data/CCGG_ClassificationSource_GeneSummary.tsv
 /Users/jvelde/github/maven/molgenis-data-cadd/data/GoldStandard_MutationTaster_TestSet.vcf
 ccgg
*
 /Users/jvelde/github/maven/molgenis-data-cadd/data/CCGG_ClassificationSource_GeneSummary.tsv
 /Users/jvelde/github/maven/molgenis-data-cadd/data/GoldStandard_UMCG_MVLs_noClinVar.vcf
 provean
*
 /Users/jvelde/github/maven/molgenis-data-cadd/data/CCGG_ClassificationSource_GeneSummary.tsv
 /Users/jvelde/github/maven/molgenis-data-cadd/data/GoldStandard_VariBench_ToleranceDS7_all.vcf.gz
 ccgg
* 
 /Users/jvelde/github/maven/molgenis-data-cadd/data/CCGG_ClassificationSource_GeneSummary.tsv
 /Users/jvelde/github/maven/molgenis-data-cadd/data/GoldStandard_UMCG_MVLs_noClinVar.vcf
 ponp2
* 
* MS Windows paths:
 C:\Users\Joeri\github\molgenis-data-cadd\data\CCGG_ClassificationSource_GeneSummary.tsv
 C:\Users\Joeri\github\molgenis-data-cadd\data\GoldStandard_UMCG_MVLs_noClinVar.vcf
 ccgg
* 
* @author jvelde
*
*/
public class Step9_Validation
{
	public static void main(String[] args) throws Exception
	{
		new Step9_Validation(args[0], args[1], args[2]);
	}
	
	CCGGUtils ccgg;
	HashMap<String, List<JudgedVariant>> judgedMVLVariants = new HashMap<String, List<JudgedVariant>>();

	/**
	 * Take MVL (annotated with CADD, ExAC, SnpEff)
	 * Take CCGG thresholds
	 * Check if classification matches
	 * @throws Exception 
	 */
	public Step9_Validation(String ccggLoc, String mvlLoc, String mode) throws Exception
	{
		List<String> modes = Arrays.asList(new String[]{"OurTool", "OurToolGenomeWide", "PONP2", "CADD", "MutationTaster2", "PROVEAN", "SIFT", "PolyPhen2", "MSC", "Condel"});
		if(!modes.contains(mode))
		{
			throw new Exception("mode needs to be one of : " + modes.toString());
		}
		File mvlFile = new File(mvlLoc);
		if(!mvlFile.exists())
		{
			throw new Exception("MVL file "+mvlFile+" does not exist or is directory");
		}
		loadCCGG(ccggLoc);
		scanMVL(mvlFile, mode);
		ProcessJudgedVariantMVLResults.printResults(judgedMVLVariants, mode, mvlFile.getName());
	}
	
	public void scanMVL(File mvlFile, String mode) throws Exception
	{
		
		VcfRepository vcfRepo = new VcfRepository(mvlFile, "mvl");
		
		java.util.Iterator<Entity> vcfRepoIter = vcfRepo.iterator();
		
		/**
		 * Here, we let PON-P2 or MutationTaster2 handle the classification and see what happens
		 * Obviously, before this works, you must score your variant list with PON-P2 or MutationTaster2
		 * And then link this file below
		 */
		PONP2Results p2r = null;
		MutationTaster2Results m2r = null;
		ProveanAndSiftResults ps2r = null;
		PolyPhen2Results pf2r = null;
		MSCResults mscr = null;
		CondelResults condelr = null;
		if (mode.equals("PONP2"))
		{
			p2r = new PONP2Results(new File("/Users/jvelde/Desktop/clinvarcadd/combined_datasets_for_external_scoring/ponp2/5832819162656_prediction.txt"));
		}
		if (mode.equals("MutationTaster2"))
		{
			m2r = new MutationTaster2Results(new File("/Users/jvelde/github/maven/molgenis-data-cadd/data/Prediction_UMCG_MVLs_noClinVar_MutationTaster2_output_minimized_conflictsRemoved.tsv"));
		}
		if (mode.equals("PROVEAN") || mode.equals("SIFT"))
		{
			ps2r = new ProveanAndSiftResults(new File("/Users/jvelde/Desktop/clinvarcadd/combined_datasets_for_external_scoring/provean/PREDICTIONS_minimized.txt"));
		}
		if (mode.equals("PolyPhen2"))
		{
			pf2r = new PolyPhen2Results(new File("/Users/jvelde/Desktop/clinvarcadd/combined_datasets_for_external_scoring/polyphen/POLYPHENRESULT-clean.txt"));
		}
		if (mode.equals("MSC"))
		{
			mscr = new MSCResults(new File("/Users/jvelde/github/maven/molgenis-data-cadd/data/MSC_CADD_cutoffs_ClinVar95CI.tsv"));
		}
		if (mode.equals("Condel"))
		{
			condelr = new CondelResults(new File("/Users/jvelde/Desktop/clinvarcadd/combined_datasets_for_external_scoring/condel/condel-results-cleaned.tsv"));
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
			
			Double getMAF = CCGGUtils.getInfoForAllele(record, "EXAC_AF", alt);
			double MAF = getMAF == null ? 0 : getMAF;
			Double CADDscore = CCGGUtils.getInfoForAllele(record, "CADD_SCALED", alt);
			String ann = record.getString("ANN");
			Set<String> genes = CCGGUtils.getGenesFromAnn(ann);
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
					Impact impact = CCGGUtils.getImpact(ann, gene, alt);

					Judgment judgment;
					if(mode.equals("OurToolGenomeWide"))
					{
						judgment = ccgg.naiveClassifyVariant(gene, MAF, impact, CADDscore);
					}
					else if (mode.equals("PONP2"))
					{
						judgment = p2r.classifyVariantUsingPONP2Results(chr, pos, ref, alt);
					}
					else if (mode.equals("MutationTaster2"))
					{
						judgment = m2r.classifyVariantUsingMutationTaster2Results(chr, pos, ref, alt);
					}
					else if (mode.equals("PROVEAN"))
					{
						judgment = ps2r.classifyVariantUsingProveanResults(chr, pos, ref, alt);
					}
					else if (mode.equals("SIFT"))
					{
						judgment = ps2r.classifyVariantUsingSiftResults(chr, pos, ref, alt);
					}
					else if (mode.equals("PolyPhen2"))
					{
						judgment = pf2r.classifyVariantUsingPolyPhen2Results(chr, pos, ref, alt);
					}
					else if (mode.equals("MSC"))
					{
						judgment = mscr.classifyVariantUsingMSCResults(gene, CADDscore);
					}
					else if (mode.equals("Condel"))
					{
						judgment = condelr.classifyVariantUsingCondelResults(chr, pos, ref, alt);
					}
					else if (mode.equals("CADD"))
					{
						if(CADDscore != null && CADDscore > 15)
						{
							judgment = new Judgment(Classification.Pathogn, Method.calibrated, "CADD score > 15");
						}
						else if(CADDscore != null && CADDscore <= 15)
						{
							judgment = new Judgment(Classification.Benign, Method.calibrated, "CADD score <= 15");
						}
						else{
							judgment = new Judgment(Classification.VOUS, Method.calibrated, "CADD score not available");
						}
					}
					else
					{
						//"ccgg" mode, the default
						judgment = ccgg.classifyVariant(gene, MAF, impact, CADDscore);
					}
					multipleJudgments.add(judgment);
					//addToMVLResults(judgment, mvlClassfc, mvlName, record);
				}
			}
			if(hasGeneId && !geneToIdMatchFound)
			{
				System.out.println("WARNING: bad data for variant " + chr + ":" + pos + " " + ref + "/" + alt + ", no match from ID field gene to snpeff annotations!");
				multipleJudgments.add(new Judgment(Classification.VOUS, Method.calibrated, "Bad data!"));
				//throw new Exception("no gene to id match for " + record.toString());
			}
			
			//if no judgment, add null for this variant
			if(multipleJudgments.size() == 0)
			{
				throw new Exception("No judgments! should not occur.");
			//	addToMVLResults(null, mvlClassfc, mvlName, record);
		//		continue;
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
				if(judgment.getClassification().equals(Classification.Pathogn))
				{
					nrOfPathognClsf++;
				}
				if(judgment.getConfidence().equals(Method.calibrated))
				{
					hasCalibratedJudgment = true;
				}
			}
			
			//check if we have any conflicts
			if(nrOfBenignClsf > 0 && nrOfPathognClsf > 0)
			{
				System.out.println("WARNING: conflicting classification! adding no judgment for this variant: " + chr + ":" + pos + " " + ref + "/" + alt + ", judged: " + multipleJudgments.toString() );
				addToMVLResults(new Judgment(Classification.VOUS, Method.calibrated, "Conflict!"), mvlClassfc, mvlName, record);
			}
			else
			{
				for(Judgment judgment : multipleJudgments)
				{
					//if we know we have calibrated results, wait for it, then add it, and then break
					if(hasCalibratedJudgment && judgment.getConfidence().equals(Method.calibrated))
					{
						addToMVLResults(judgment, mvlClassfc, mvlName, record);
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
	
	public void loadCCGG(String ccggLoc) throws Exception
	{
		File ccggFile = new File(ccggLoc);
		if(!ccggFile.exists())
		{
			throw new Exception("CCGG file "+ccggFile+" does not exist or is directory");
		}
		CCGGUtils ccggUtils = new CCGGUtils(ccggFile);
		this.ccgg = ccggUtils;
	}
}
