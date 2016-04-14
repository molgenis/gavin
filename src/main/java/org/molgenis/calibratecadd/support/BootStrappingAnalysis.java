package org.molgenis.calibratecadd.support;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Set;

import org.molgenis.calibratecadd.Step9_Validation;
import org.molgenis.calibratecadd.Step9_Validation.ToolNames;
import org.molgenis.data.Entity;
import org.molgenis.data.annotation.entity.impl.SnpEffAnnotator.Impact;
import org.molgenis.data.annotation.joeri282exomes.CCGGUtils;
import org.molgenis.data.annotation.joeri282exomes.Judgment;
import org.molgenis.data.annotation.joeri282exomes.Judgment.Classification;
import org.molgenis.data.annotation.joeri282exomes.Judgment.Method;
import org.molgenis.data.vcf.VcfRepository;

public class BootStrappingAnalysis
{

	public static void main(String[] args) throws Exception
	{

		new BootStrappingAnalysis(args[0], args[1], Integer.parseInt(args[2]));

	}
	
	int TP = 0;
	int TN = 0;
	int FP = 0;
	int FN = 0;
	int VOUS = 0;
	int judgmentsInCalibratedGenes = 0;

	/**
	 * checked to give the exact same behaviour as Step9_Validation
	 * Except now we do random subsets of variants
	 * @param vcfFile
	 * @param gavinFile
	 * @throws Exception
	 */
	public BootStrappingAnalysis(String vcfFile, String gavinFile, int bootstrapsize) throws Exception
	{

	
		CCGGUtils gavin = Step9_Validation.loadCCGG(gavinFile);
		
		// file with combined variants has 25,995 variants
		File variantList = new File(vcfFile);

		// get 100 random numbers in the range 1-25995
		ArrayList<Integer> list = new ArrayList<Integer>();
		ArrayList<Integer> randomX = new ArrayList<Integer>();
		for (int i = 1; i <= 25995; i++)
		{
			list.add(i);
		}
		Collections.shuffle(list);
		for (int i = 0; i < bootstrapsize; i++)
		{
			randomX.add(list.get(i));
		}
		
		System.out.println("selected line nrs: " + randomX.toString());

		VcfRepository vcfRepo = new VcfRepository(variantList, "vcf");

		java.util.Iterator<Entity> vcfRepoIter = vcfRepo.iterator();
		
		int lineNr = 0;
		while (vcfRepoIter.hasNext())
		{
			lineNr++;
			
			Entity record = vcfRepoIter.next();
		
			if(!randomX.contains(lineNr))
			{
				continue;
			}
			System.out.println("get this line: " + lineNr);
			
			
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
					Judgment judgment = gavin.classifyVariant(gene, MAF, impact, CADDscore);
					multipleJudgments.add(judgment);
				}
			}
			if(hasGeneId && !geneToIdMatchFound)
			{
	//			System.out.println("WARNING: bad data for variant " + chr + ":" + pos + " " + ref + "/" + alt + ", no match from ID field gene to snpeff annotations!");
				multipleJudgments.add(new Judgment(Classification.VOUS, Method.calibrated, "Bad data!"));
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
				if(judgment.getClassification().equals(Classification.Pathogn))
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
				count(Classification.VOUS, mvlClassfc);
			//	System.out.println("WARNING: conflicting classification! adding no judgment for this variant: " + chr + ":" + pos + " " + ref + "/" + alt + ", judged: " + multipleJudgments.toString() );
			}
			else
			{
				for(Judgment judgment : multipleJudgments)
				{
					//if we know we have calibrated results, wait for it, then add it, and then break
					if(hasCalibratedJudgment && judgment.getConfidence().equals(Method.calibrated))
					{
					//	addToMVLResults(judgment, mvlClassfc, mvlName, record);
						judgmentsInCalibratedGenes++;
						count(judgment.getClassification(), mvlClassfc);
						break;
					}
					//if not, might as well add this one and be done
					//TODO: this means there may be multiple verdicts, e.g. 2x BENIGN for context in two genes, but we only add 1 of them, to keep things a bit more simple
					else if(!hasCalibratedJudgment)
					{
					//	addToMVLResults(judgment, mvlClassfc, mvlName, record);
						
						count(judgment.getClassification(), mvlClassfc);
						break;
					}
				}
			}

		}
		
		
		System.out.println("TN = " + TN);
		System.out.println("TP = " + TP);
		System.out.println("FP = " + FP);
		System.out.println("FN = " + FN);
		System.out.println("VOUS = " + VOUS);
		double MCC = ProcessJudgedVariantMVLResults.getMCC(TP, TN, FP, FN);
		System.out.println("MCC = " + MCC);
		double cov = (TP+TN+FP+FN)/(double)(TP+TN+FP+FN+VOUS);
		System.out.println("MCCcovadj = " + (cov*MCC));
		double percCalib = judgmentsInCalibratedGenes/(double)bootstrapsize;
		System.out.println("% calibrated: " + judgmentsInCalibratedGenes + "/" + bootstrapsize + " = " + percCalib);
		
	}
	
	public void count (Classification observed, String expected)
	{
		if(observed.equals(Classification.Benign) && (expected.equals("B") ||  expected.equals("LB")))
		{
			TN++;
		}
		
		if(observed.equals(Classification.Benign) && (expected.equals("P") ||  expected.equals("LP")))
		{
			FN++;
		}
		
		if(observed.equals(Classification.Pathogn) && (expected.equals("B") ||  expected.equals("LB")))
		{
			FP++;
		}
		
		if(observed.equals(Classification.Pathogn) && (expected.equals("P") ||  expected.equals("LP")))
		{
			TP++;
		}
		
		if(observed.equals(Classification.VOUS) && (expected.equals("P") ||  expected.equals("LP") || expected.equals("B") ||  expected.equals("LB")))
		{
			VOUS++;
		}
				
		
	}

}
