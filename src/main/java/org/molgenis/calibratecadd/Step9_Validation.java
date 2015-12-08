package org.molgenis.calibratecadd;

import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Set;

import org.molgenis.calibratecadd.support.PONP2Results;
import org.molgenis.calibratecadd.support.VariantClassificationException;
import org.molgenis.calibratecadd.support.JudgedVariant;
import org.molgenis.calibratecadd.support.JudgedVariant.ExpertClassification;
import org.molgenis.calibratecadd.support.ProcessJudgedVariantMVLResults;
import org.molgenis.data.Entity;
import org.molgenis.data.annotation.entity.impl.SnpEffAnnotator.Impact;
import org.molgenis.data.annotation.joeri282exomes.CCGGUtils;
import org.molgenis.data.annotation.joeri282exomes.Judgment;
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
		if(!mode.equals("full") && !mode.equals("naiveonly") && !mode.equals("ponp2"))
		{
			throw new Exception("mode needs to be 'full', 'naiveonly', or 'ponp2'");
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
			String geneFromId = id.split(":", -1)[0];
			String mvlClassfc = record.getString("CLSF");
			String mvlName = record.getString("MVL");
			
			boolean geneToIdMatch = false;
			for(String gene : genes)
			{
				if(gene.equals(geneFromId))
				{
					geneToIdMatch = true;
					Impact impact = CCGGUtils.getImpact(ann, gene, alt);
					try
					{
						Judgment judgment;
						if(mode.equals("naiveonly"))
						{
							judgment = ccgg.naiveClassifyVariant(gene, MAF, impact, CADDscore);
						}
						else if (mode.equals("ponp2"))
						{
							PONP2Results p2r = new PONP2Results(new File("/Users/jvelde/github/maven/molgenis-data-cadd/data/PONP2_predictions.txt"));
							judgment = p2r.classifyVariantUsingPONP2Results(chr, pos, ref, alt);
						}
						else
						{
							judgment = ccgg.classifyVariant(gene, MAF, impact, CADDscore);
						}
						addToMVLResults(judgment, mvlClassfc, mvlName, record);
					}
					catch(VariantClassificationException e)
					{
						addToMVLResults(null, mvlClassfc, mvlName, record);
					}
					break;
				}
			}
			if(!geneToIdMatch)
			{
				throw new Exception("no gene to id match for " + record.toString());
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
