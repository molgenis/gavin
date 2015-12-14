package org.molgenis.calibratecadd;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Set;

import org.molgenis.calibratecadd.support.MutationTaster2Results;
import org.molgenis.calibratecadd.support.PONP2Results;
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
		List<String> modes = Arrays.asList(new String[]{"ccgg", "naive", "ponp2", "mutationtaster2", "provean", "sift"});
		if(!modes.contains(mode))
		{
			throw new Exception("mode needs to be one of : " + modes.toString());
		}
		loadCCGG(ccggLoc);
		scanMVL(mvlLoc, mode);
		ProcessJudgedVariantMVLResults.printResults(judgedMVLVariants);
	}
	
	public void scanMVL(String mvlLoc, String mode) throws Exception
	{
		File mvlFile = new File(mvlLoc);
		if(!mvlFile.exists())
		{
			throw new Exception("MVL file "+mvlFile+" does not exist or is directory");
		}
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
		if (mode.equals("ponp2"))
		{
			p2r = new PONP2Results(new File("/Users/jvelde/github/maven/molgenis-data-cadd/data/PONP2_predictions.txt"));
		}
		if (mode.equals("mutationtaster2"))
		{
			m2r = new MutationTaster2Results(new File("/Users/jvelde/github/maven/molgenis-data-cadd/data/MutationTaster2_output_minimized.tsv"));
		}
		if (mode.equals("provean") || mode.equals("sift"))
		{
			ps2r = new ProveanAndSiftResults(new File("/Users/jvelde/github/maven/molgenis-data-cadd/data/provean_sift_predictions_minimized.tsv"));
		}

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
					try
					{
						Judgment judgment;
						if(mode.equals("naive"))
						{
							judgment = ccgg.naiveClassifyVariant(gene, MAF, impact, CADDscore);
						}
						else if (mode.equals("ponp2"))
						{
							judgment = p2r.classifyVariantUsingPONP2Results(chr, pos, ref, alt);
						}
						else if (mode.equals("mutationtaster2"))
						{
							judgment = m2r.classifyVariantUsingMutationTaster2Results(chr, pos, ref, alt);
						}
						else if (mode.equals("provean"))
						{
							judgment = ps2r.classifyVariantUsingProveanResults(chr, pos, ref, alt);
						}
						else if (mode.equals("sift"))
						{
							judgment = ps2r.classifyVariantUsingSiftResults(chr, pos, ref, alt);
						}
						else
						{
							//"ccgg" mode, the default
							judgment = ccgg.classifyVariant(gene, MAF, impact, CADDscore);
						}
						multipleJudgments.add(judgment);
						//addToMVLResults(judgment, mvlClassfc, mvlName, record);
					}
					catch(VariantClassificationException e)
					{
						//addToMVLResults(null, mvlClassfc, mvlName, record);
					}
				}
			}
			if(hasGeneId && !geneToIdMatchFound)
			{
				System.out.println("WARNING: bad data for variant " + chr + ":" + pos + " " + ref + "/" + alt + ", no match from ID field gene to snpeff annotations!");
				//throw new Exception("no gene to id match for " + record.toString());
			}
			
			//if no judgment, add null for this variant
			if(multipleJudgments.size() == 0)
			{
				addToMVLResults(null, mvlClassfc, mvlName, record);
				continue;
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
				System.out.println("WARNING: conflicting classification! not adding this variant: " + chr + ":" + pos + " " + ref + "/" + alt + ", judged: " + multipleJudgments.toString() );
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
