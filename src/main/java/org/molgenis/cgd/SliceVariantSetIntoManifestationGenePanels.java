package org.molgenis.cgd;

import net.didion.jwnl.data.Exc;
import org.molgenis.calibratecadd.support.GavinUtils;
import org.molgenis.data.Entity;
import org.molgenis.data.vcf.VcfRepository;
import org.molgenis.data.vcf.utils.VcfWriterUtils;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.*;

/**
 * Created by joeri on 6/2/16.
 */
public class SliceVariantSetIntoManifestationGenePanels {

    public SliceVariantSetIntoManifestationGenePanels(File vcfFile, File cgdFile, File outputDir) throws Exception {

        VcfRepository vcf = new VcfRepository(vcfFile, "vcf");
        HashMap<String, CGDEntry> cgd = LoadCGD.loadCGD(cgdFile);
        HashMap<String, List<Entity>> manifestionToRecord = new HashMap<>();
        HashMap<String, Set<String>> manifestionToGenes = new HashMap<>();

        if(!outputDir.isDirectory())
        {
            throw new Exception("output location is not a directory");
        }

        Iterator<Entity> it = vcf.iterator();

        while (it.hasNext()) {
            Entity record = it.next();
            String ann = record.getString("ANN");
            Set<String> genes = GavinUtils.getGenesFromAnn(ann);

            Set<String> manifestationCategoriesAlreadyWritten = new HashSet<>();

            boolean inCGD = false;
            for(String gene : genes)
            {
                if(cgd.keySet().contains(gene))
                {
                    inCGD = true;
                    if(cgd.get(gene).getManifestationCategoriesList().size() == 0)
                    {
                        throw new Exception("gene in CGD but no ManifestationCategories: " + gene );
                    }
                    for(String manifestCat : cgd.get(gene).getManifestationCategoriesList())
                    {
                        if(manifestationCategoriesAlreadyWritten.contains(manifestCat))
                        {
                            //prevent writing out to the same manifestation category more than once when there are multiple CGD genes that map to the same
                            System.out.println("already written variant to manifestation category " + manifestCat + ", skipping!");
                            continue;
                        }
                        manifestationCategoriesAlreadyWritten.add(manifestCat);

                        List<Entity> variantList = manifestionToRecord.get(manifestCat);
                        if(variantList == null) {
                            variantList = new ArrayList<>();
                            manifestionToRecord.put(manifestCat, variantList);
                        }
                        variantList.add(record);

                        Set<String> geneList = manifestionToGenes.get(manifestCat);
                        if(geneList == null) {
                            geneList = new HashSet<>();
                            manifestionToGenes.put(manifestCat, geneList);
                        }
                        geneList.add(gene);

                    }
                }
            }

            if(!inCGD)
            {
                //only when 0 of the genes is found in CGD we consider a variant to be 'not in CGD'
                // this prevents the sitation  where e.g. TTN variants that also lie in MIR548N and TTN-AS1 pseudogene go into "Cardiovascular" as well as "NotInCGD" ... which is a bit weird
                List<Entity> variantList = manifestionToRecord.get("NotInCGD");
                if(variantList == null) {
                    variantList = new ArrayList<>();
                    manifestionToRecord.put("NotInCGD", variantList);
                }
                variantList.add(record);

                Set<String> geneList = manifestionToGenes.get("NotInCGD");
                if(geneList == null) {
                    geneList = new HashSet<>();
                    manifestionToGenes.put("NotInCGD", geneList);
                }
                geneList.addAll(genes);
            }

        }

        for(String key : manifestionToRecord.keySet())
        {
            System.out.println(key + " has " + manifestionToRecord.get(key).size() + " variants in " + manifestionToGenes.get(key).size() + " genes");
            writeToVcf(manifestionToRecord.get(key), new File(outputDir, key.replace("/", "_") + "")); //dont add VCF for convenience but maybe you want to
        }

    }

    public void writeToVcf(List<Entity> variants, File writeTo) throws IOException {
        String header = "##fileformat=VCFv4.0\n" +
                "##INFO=<ID=CADD,Number=.,Type=Float,Description=\"na\">\n" +
                "##INFO=<ID=CADD_SCALED,Number=.,Type=Float,Description=\"na\">\n" +
                "##VEP=v82 cache=/data_ensembl/vep/grch37/release-83/homo_sapiens/83_GRCh37 db=homo_sapiens_core_83_37@ensdb-web-14 sift=sift5.2.2 polyphen=2.2.2 COSMIC=71 ESP=20141103 gencode=GENCODE 19 HGMD-PUBLIC=20152 genebuild=2011-04 regbuild=13 ClinVar=201507 dbSNP=144 assembly=GRCh37.p13\n" +
                "##INFO=<ID=MVL,Number=.,Type=String,Description=\"The MVL this variant belongs to\">\n" +
                "##INFO=<ID=CLSF,Number=.,Type=String,Description=\"The classification of this variant\">\n" +
                "##INFO=<ID=CSQ,Number=.,Type=String,Description=\"Consequence annotations from Ensembl VEP. Format: Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|DISTANCE|STRAND|SYMBOL_SOURCE|HGNC_ID|TSL|APPRIS|SIFT|PolyPhen|GMAF|AFR_MAF|AMR_MAF|EAS_MAF|EUR_MAF|SAS_MAF|AA_MAF|EA_MAF|ExAC_MAF|ExAC_Adj_MAF|ExAC_AFR_MAF|ExAC_AMR_MAF|ExAC_EAS_MAF|ExAC_FIN_MAF|ExAC_NFE_MAF|ExAC_OTH_MAF|ExAC_SAS_MAF|CLIN_SIG|SOMATIC|PHENO|PUBMED|MOTIF_NAME|MOTIF_POS|HIGH_INF_POS|MOTIF_SCORE_CHANGE\">\n" +
                "##SnpEffVersion=\"4.2 (build 2015-12-05), by Pablo Cingolani\"\n" +
                "##SnpEffCmd=\"SnpEff  hg19 -noStats -lof ../onco/oncovariants.fix.vcf \"\n" +
                "##INFO=<ID=ANN,Number=.,Type=String,Description=\"Functional annotations: 'Allele | Annotation | Annotation_Impact | Gene_Name | Gene_ID | Feature_Type | Feature_ID | Transcript_BioType | Rank | HGVS.c | HGVS.p | cDNA.pos / cDNA.length | CDS.pos / CDS.length | AA.pos / AA.length | Distance | ERRORS / WARNINGS / INFO' \">\n" +
                "##INFO=<ID=LOF,Number=.,Type=String,Description=\"Predicted loss of function effects for this variant. Format: 'Gene_Name | Gene_ID | Number_of_transcripts_in_gene | Percent_of_transcripts_affected' \">\n" +
                "##INFO=<ID=NMD,Number=.,Type=String,Description=\"Predicted nonsense mediated decay effects for this variant. Format: 'Gene_Name | Gene_ID | Number_of_transcripts_in_gene | Percent_of_transcripts_affected' \">\n" +
                "##INFO=<ID=EXAC_AF,Number=.,Type=String,Description=\"The ExAC allele frequency\">\n" +
                "##INFO=<ID=EXAC_AC_HOM,Number=.,Type=String,Description=\"The ExAC homozygous alternative genotype count\">\n" +
                "##INFO=<ID=EXAC_AC_HET,Number=.,Type=String,Description=\"The ExAC heterozygous genotype count\">\n" +
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";

        FileWriter fw = new FileWriter(writeTo);
        BufferedWriter bw = new BufferedWriter(fw);
        bw.write(header);
        for(Entity e : variants)
        {
            VcfWriterUtils.writeToVcf(e, bw);
            bw.write("\n");
        }
        bw.close();
    }

    public static void main (String[] args) throws Exception {

        if(args.length != 3)
        {
            throw new Exception(("please provide: VCF location, CGD location, output directory"));
        }
       new SliceVariantSetIntoManifestationGenePanels(new File(args[0]), new File(args[1]), new File(args[2]));

    }

}
