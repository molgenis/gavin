package org.molgenis.cgd;

import java.util.ArrayList;
import java.util.List;

public class CGDEntry
{
	private String gene;
	private String hgnc_id;
	private String entrez_gene_id;
	private String condition;
	private String inheritance;
	private String age_group;
	private String allelicConditions;
	private String manifestationCategories;
	private List<String> manifestationCategoriesList;
	private String interventionCategories;
	private String comments;
	private String interventionOrRationale;
	private String references;
	private generalizedInheritance generalizedInheritance; //used for gene inheritance matching!


	public enum generalizedInheritance {

		DOMINANT, RECESSIVE, DOMINANT_OR_RECESSIVE, XL_LINKED, BLOODGROUP, OTHER, NOTINCGD;

		public static boolean hasKnownInheritance(generalizedInheritance gi) {
			return (gi == RECESSIVE || gi == DOMINANT || gi == DOMINANT_OR_RECESSIVE || gi == XL_LINKED) ? true : false;
		}
	}
	
	public CGDEntry(String gene, String hgnc_id, String entrez_gene_id,
			String condition, String inheritance, String age_group,
			String allelicConditions, String manifestationCategories,
			String interventionCategories, String comments,
			String interventionOrRationale, String references) {
		super();
		this.gene = gene;
		this.hgnc_id = hgnc_id;
		this.entrez_gene_id = entrez_gene_id;
		this.condition = condition;
		this.inheritance = inheritance;
		this.age_group = age_group;
		this.allelicConditions = allelicConditions;
		this.manifestationCategories = manifestationCategories;
		this.interventionCategories = interventionCategories;
		this.comments = comments;
		this.interventionOrRationale = interventionOrRationale;
		this.references = references;
		ArrayList<String> manifestationCategoriesList = new ArrayList<>();
		for(String manifestCat: manifestationCategories.split(";", -1))
		{
			if(!manifestCat.trim().equals("")) {
				manifestationCategoriesList.add(manifestCat.trim());
			}
		}
		this.manifestationCategoriesList = manifestationCategoriesList;
	}

	public List<String> getManifestationCategoriesList() {
		return manifestationCategoriesList;
	}

	public String getGene() {
		return gene;
	}



	public String getHgnc_id() {
		return hgnc_id;
	}



	public String getEntrez_gene_id() {
		return entrez_gene_id;
	}



	public String getCondition() {
		return condition;
	}



	public String getInheritance() {
		return inheritance;
	}



	public String getAge_group() {
		return age_group;
	}



	public String getAllelicConditions() {
		return allelicConditions;
	}



	public String getManifestationCategories() {
		return manifestationCategories;
	}



	public String getInterventionCategories() {
		return interventionCategories;
	}



	public String getComments() {
		return comments;
	}



	public String getInterventionOrRationale() {
		return interventionOrRationale;
	}



	public String getReferences() {
		return references;
	}



	public generalizedInheritance getGeneralizedInheritance() {
		return generalizedInheritance;
	}



	public void setGeneralizedInheritance(generalizedInheritance generalizedInheritance) {
		this.generalizedInheritance = generalizedInheritance;
	}



	@Override
	public String toString() {
		return gene + "\t" + hgnc_id
				+ "\t" + entrez_gene_id + "\t"
				+ condition + "\t" + inheritance + "\t"
				+ age_group + "\t" + allelicConditions
				+ "\t" + manifestationCategories
				+ "\t" + interventionCategories
				+ "\t" + comments + "\t"
				+ interventionOrRationale + "\t" + references
				+ "\t" + generalizedInheritance;
	}
	
	
	
}
