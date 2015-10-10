package org.molgenis.calibratecadd;

public class Step5_CADDScores 
{
	/**
	 * Run the output file of the previous step through the CADD web service at
	 * cadd.gs.washington.edu/score
	 * 
	 * The filenames ends in '.cadd', e.g. clinvar.patho.fix.snpeff.exac.vcf.cadd
	 * 
	 * You may have to split up the file in smaller ones, say aroung 50k variants each or 1 MB in size.
	 * 
	 * When done, combine the result files (for download at the URLs provided, after a certain waiting period)
	 * back into 1 file and proceed to the next step.
	 * 
	 * 
	 * Alternatively, ofcourse, run CADD locally. But keep in mind that you'll have to run against multiple CADD
	 * files to get all SNVs and indels scores as best as possible. Ideally, install the CADD algorithm for live
	 * (and complete) calculation of scores.
	 * 
	 */
}
