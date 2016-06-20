package org.molgenis.calibratecadd.support;

import org.molgenis.data.Entity;
import org.molgenis.data.annotation.entity.impl.gavin.Judgment;
import org.molgenis.data.annotation.entity.impl.gavin.Judgment.Classification;
import org.molgenis.data.annotation.entity.impl.gavin.Judgment.Method;
import org.molgenis.data.vcf.VcfRepository;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.HashMap;
import java.util.Scanner;

public class PredictSNP2Results
{

	//chr:pos,ref,alt to PredictSNP2 consensus ("PSNPE" info field)
	HashMap<String, HashMap<String, String>> toolToEstimatedEffect;

	/**
	 *
	 */
	public PredictSNP2Results(File ps2Results) throws IOException
	{
		VcfRepository vcfRepo = new VcfRepository(ps2Results, "mvl");
		HashMap<String, String> PredictSNP2 = new HashMap<String, String>();
		HashMap<String, String> FATHMM = new HashMap<String, String>();
		HashMap<String, String> GWAVA = new HashMap<String, String>();
		HashMap<String, String> FunSeq = new HashMap<String, String>();
		HashMap<String, String> DANN = new HashMap<String, String>();
		toolToEstimatedEffect = new HashMap<String, HashMap<String, String>>();

		java.util.Iterator<Entity> vcfRepoIter = vcfRepo.iterator();

		while (vcfRepoIter.hasNext()) {
			Entity record = vcfRepoIter.next();

			String chr = record.getString("#CHROM");
			String pos = record.getString("POS");
			String ref = record.getString("REF");
			String alt = record.getString("ALT");
			String PSNPE = record.getString("PSNPE");
			String FATE = record.getString("FATE");
			String GWAVAE = record.getString("GWAVAE");
			String DANNE = record.getString("DANNE");
			String FUNE = record.getString("FUNE");
			PredictSNP2.put(chr + "_" + pos + "_" + ref + "_"+alt,PSNPE);
			FATHMM.put(chr + "_" + pos + "_" + ref + "_"+alt,FATE);
			GWAVA.put(chr + "_" + pos + "_" + ref + "_"+alt,GWAVAE);
			FunSeq.put(chr + "_" + pos + "_" + ref + "_"+alt,FUNE);
			DANN.put(chr + "_" + pos + "_" + ref + "_"+alt,DANNE);
		}
		toolToEstimatedEffect.put("PredictSNP2", PredictSNP2);
		toolToEstimatedEffect.put("FATHMM", FATHMM);
		toolToEstimatedEffect.put("GWAVA", GWAVA);
		toolToEstimatedEffect.put("FunSeq", FunSeq);
		toolToEstimatedEffect.put("DANN", DANN);
	}

	public static void main (String[] args) throws IOException {
		new PredictSNP2Results(new File("/Users/joeri/github/gavin/data/predictions", "PredictSNP2.vcf.gz"));
	}
	
	public Judgment classifyVariantUsingPredictSNP2Results(String chr, String pos, String ref, String alt, String tool) throws Exception
	{
		String key = chr+"_"+pos+"_"+ref+"_"+alt;
		if(toolToEstimatedEffect.get(tool).containsKey(key))
		{
			String consensus = toolToEstimatedEffect.get(tool).get(key);

			if(consensus.contains("neutral"))
			{
				return new Judgment(Classification.Benign, Method.calibrated, key, tool + " result 'neutral'");
			}
			else if(consensus.contains("deleterious"))
			{
				return new Judgment(Classification.Pathogenic, Method.calibrated, key, tool + " result 'deleterious'");
			}
			else
			{
				return new Judgment(Classification.VOUS, Method.calibrated, key, "Unknown "+tool+" result: " + consensus);
			}
		}
		else
		{
			return new Judgment(Classification.VOUS, Method.calibrated, key, "No "+tool+" result");
		}
	}
	
	
}
