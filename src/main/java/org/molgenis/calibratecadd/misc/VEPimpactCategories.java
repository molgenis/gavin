package org.molgenis.calibratecadd.misc;


public class VEPimpactCategories
{
	
	public static int getImpactRank(String impact) throws Exception
	{
		switch (impact) {
			case "HIGH":  return 3;
			case "MODERATE":  return 2;
			case "LOW":  return 1;
			case "MODIFIER":  return 0;
			default: throw new Exception("impact unknown: " + impact);
		}
	}

	//TODO these definitions seem slightly different than SnpEffs ?
	//check and harmonize with http://snpeff.sourceforge.net/SnpEff_manual.html !
	public static String getImpact(String consequence) throws Exception
	{
		switch (consequence) {
			case "transcript_ablation":  return "HIGH";
			case "splice_acceptor_variant":  return "HIGH";
			case "splice_donor_variant":  return "HIGH";
			case "stop_gained":  return "HIGH";
			case "frameshift_variant":  return "HIGH";
			case "stop_lost":  return "HIGH";
			case "start_lost":  return "HIGH";
			case "transcript_amplification":  return "HIGH";
			
			case "inframe_insertion":  return "MODERATE";
			case "inframe_deletion":  return "MODERATE";
			case "missense_variant":  return "MODERATE";
			case "protein_altering_variant":  return "MODERATE";
			
			case "splice_region_variant":  return "LOW";
			case "incomplete_terminal_codon_variant":  return "LOW";
			case "stop_retained_variant":  return "LOW";
			case "synonymous_variant":  return "LOW";
			
			case "coding_sequence_variant":  return "MODIFIER";
			case "mature_miRNA_variant":  return "MODIFIER";
			case "5_prime_UTR_variant":  return "MODIFIER";
			case "non_coding_transcript_exon_variant":  return "MODIFIER";
			case "intron_variant":  return "MODIFIER";
			case "NMD_transcript_variant":  return "MODIFIER";
			case "non_coding_transcript_variant":  return "MODIFIER";
			case "upstream_gene_variant":  return "MODIFIER";
			case "downstream_gene_variant":  return "MODIFIER";
			case "TFBS_ablation":  return "MODIFIER";
			case "TFBS_amplification":  return "MODIFIER";
			case "TF_binding_site_variant":  return "MODIFIER";
			case "regulatory_region_ablation":  return "MODIFIER";
			case "regulatory_region_amplification":  return "MODIFIER";
			case "feature_elongation":  return "MODIFIER";
			case "regulatory_region_variant":  return "MODIFIER";
			case "feature_truncation":  return "MODIFIER";
			case "intergenic_variant":  return "MODIFIER";
			
            default: throw new Exception("consequence unknown: " + consequence);
		}
                 
	}
	
	
	
}
