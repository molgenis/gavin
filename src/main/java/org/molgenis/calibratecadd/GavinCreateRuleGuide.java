package org.molgenis.calibratecadd;

import org.molgenis.calibratecadd.support.GavinUtils;
import org.molgenis.data.annotation.entity.impl.gavin.GavinAlgorithm;
import org.molgenis.data.annotation.entity.impl.gavin.GavinEntry;
import org.molgenis.data.annotation.entity.impl.gavin.GavinEntry.Category;

import java.io.File;
import java.io.PrintWriter;
import java.util.HashMap;

/**
 * Created by joeri on 9/22/16.
 */
public class GavinCreateRuleGuide {

    public static void main(String[] args) throws Exception
    {
        if(args.length != 2)
        {
            throw new Exception("please provide: gavin file (e.g. /github/gavin/data/predictions), and version, (e.g. r0.2)");
        }

        String version = args[1];
        String gavinCalibrationFile = args[0] + File.separatorChar + "GAVIN_calibrations_"+version+".tsv";
        String output = args[0] + File.separatorChar + "GAVIN_ruleguide_"+version+".tsv";
        File outputFile = new File(output);

        //dev
        outputFile.delete();

        outputFile.createNewFile();
        new GavinCreateRuleGuide(gavinCalibrationFile, outputFile).createRuleGuide(version);
    }

    HashMap<String, GavinEntry> gavinData;
    PrintWriter pw;

    public GavinCreateRuleGuide(String gavinCalibrationFile, File outputFile) throws Exception
    {
        this.gavinData = Benchmark.loadGAVIN(gavinCalibrationFile).getGeneToEntry();
        this.pw = new PrintWriter(outputFile);
    }

    public void createRuleGuide(String version) throws Exception {

        String sep = "\t";
        String NA = "n/a";

        pw.println("################################");
        pw.println("## GAVIN Applied Rule Guide " + version);
        pw.println("################################");
        pw.println("## ");
        pw.println("## Variant can be interpreted by following the columns from left to right.");
        pw.println("## This classification scheme was implemented in Java, and used to");
        pw.println("## benchmark GAVIN in the paper (Van der Velde et al., using r0.2), see:");
        pw.println("## https://github.com/molgenis/molgenis and https://github.com/molgenis/gavin");
        pw.println("## ");
        pw.println("## Explanation of the gene calibration categories:");
        pw.println("## C1 = CADD scores highly predictive (pval < 0.01)");
        pw.println("## C2 = CADD scores predictive (pval < 0.05)");
        pw.println("## C3 = CADD may be predictive (pval > 0.05 but with few samples)");
        pw.println("## C4 = CADD not so predictive (pval > 0.05 with enough samples)");
        pw.println("## C5 = CADD not so predictive (population CADD > pathogenic CADD)");
        pw.println("## I1 = HIGH impact unique for pathogenic variants");
        pw.println("## I2 = MODERATE or HIGH impact unique for pathogenic variants");
        pw.println("## I3 = LOW or MODERATE or HIGH impact unique for pathogenic variants");
        pw.println("## T1 = Too few exac variants after filtering with pathogenic 95th percentile MAF.");
        pw.println("## T2 = Too few exac variants after filtering with impact distribution.");
        pw.println("## N1 = Too few ClinVar variants for calibration at this time.");
        pw.println("## N2 = Too few ExAC variants found for calibration.");
        pw.println("## ");
        pw.println("## Genome-wide rules are applied if the gene-specific rules have failed to classify.");
        pw.println("## These general rules are as follows:");
        pw.println("## If impact equals MODIFIER -> benign");
        pw.println("## If variant MAF greater than " + GavinAlgorithm.GENOMEWIDE_MAF_THRESHOLD + " -> benign");
        pw.println("## If CADD less than "+GavinAlgorithm.GENOMEWIDE_CADD_THRESHOLD+" -> benign");
        pw.println("## If CADD greater than "+GavinAlgorithm.GENOMEWIDE_CADD_THRESHOLD+" -> pathogenic");
        pw.println("## ");
        pw.println("Gene" + sep + "Category" + sep + "PathogenicIfCADDScoreGreaterThan" + sep + "BenignIfCADDScoreLessThan" + sep + "BenignIfMAFGreaterThan" + sep + "PathogenicIfImpactEquals");

        for(String gene : gavinData.keySet())
        {
            Category category = gavinData.get(gene).category;
            Double pathoMAFThreshold = gavinData.get(gene).PathoMAFThreshold != null ?
                    gavinData.get(gene).PathoMAFThreshold * GavinAlgorithm.extraSensitivityFactor * 2 : null;
            Double meanPathogenicCADDScore = gavinData.get(gene).MeanPathogenicCADDScore != null ?
                    gavinData.get(gene).MeanPathogenicCADDScore - GavinAlgorithm.extraSensitivityFactor : null;
            Double meanPopulationCADDScore = gavinData.get(gene).MeanPopulationCADDScore != null ?
                    gavinData.get(gene).MeanPopulationCADDScore - GavinAlgorithm.extraSensitivityFactor : null;
            Double spec95thPerCADDThreshold = gavinData.get(gene).Spec95thPerCADDThreshold != null ?
                    gavinData.get(gene).Spec95thPerCADDThreshold - GavinAlgorithm.extraSensitivityFactor : null;
            Double sens95thPerCADDThreshold = gavinData.get(gene).Sens95thPerCADDThreshold != null ?
                    gavinData.get(gene).Sens95thPerCADDThreshold - GavinAlgorithm.extraSensitivityFactor : null;

            pathoMAFThreshold = pathoMAFThreshold != null ? Math.min(pathoMAFThreshold, 1) : null;

            switch (category) {
                case C1: case C2:
                    pw.println(gene + sep + category + sep + Math.max(meanPathogenicCADDScore, 0) + sep + Math.max(meanPopulationCADDScore, 0) + sep + pathoMAFThreshold + sep + NA);
                    break;
                case C3: case C4: case C5:
                    pw.println(gene + sep + category + sep + Math.max(spec95thPerCADDThreshold, 0) + sep + Math.max(sens95thPerCADDThreshold, 0) + sep + pathoMAFThreshold + sep + NA);
                    break;
                case I1: case I2: case I3:
                    pw.println(gene + sep + category + sep + NA + sep + NA + sep + pathoMAFThreshold + sep + (category == Category.I1 ? "HIGH" : category == Category.I2 ? "HIGH,MODERATE" : category == Category.I3 ? "HIGH,MODERATE,LOW" : null));
                    break;
                case T1: case T2:
                    pw.println(gene + sep + category + sep + NA + sep + NA + sep + pathoMAFThreshold + sep + NA);
                    break;
                case N1: case N2:
                    pw.println(gene + sep + category + sep + NA + sep + NA + sep + NA + sep + NA);
                    break;

                default:
                    throw new Exception("Invalid category: " + category);
            }
        }

        pw.flush();
        pw.close();
    }


}
