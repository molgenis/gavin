package org.molgenis.calibratecadd.support;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Scanner;
import java.util.Set;

import org.apache.commons.math3.stat.descriptive.moment.Mean;
import org.molgenis.data.annotation.joeri282exomes.CCGGEntry;
import org.molgenis.data.annotation.joeri282exomes.CCGGEntry.Category;
import org.molgenis.data.annotation.joeri282exomes.CCGGUtils;

public class CompareToCGDandGDI
{
	
	
	public static void main(String[] args) throws Exception
	{
		File cgdFile = new File(args[0]);
		File gdiFile = new File(args[1]);
		File ccggFile = new File(args[2]);
	
		CompareToCGDandGDI comp = new CompareToCGDandGDI(cgdFile, gdiFile, ccggFile);
		comp.loadData();
		comp.compareToCGD();
		comp.compareToGDI();
	}
	

	File cgdFile;
	File ccggFile;
	File gdiFile;
	HashMap<String, String> geneToInheritance;
	HashMap<String, String> geneToOnset;
	HashMap<String, String> geneToDamage;
	CCGGUtils cu;
	
	public CompareToCGDandGDI(File cgdFile, File gdiFile, File ccggFile) throws FileNotFoundException
	{
		if (!cgdFile.isFile())
		{
			throw new FileNotFoundException("cgdFile " + cgdFile.getAbsolutePath()
					+ " does not exist or is directory");
		}
		if (!gdiFile.isFile())
		{
			throw new FileNotFoundException("gdiFile " + gdiFile.getAbsolutePath()
					+ " does not exist or is directory");
		}
		if (!ccggFile.isFile())
		{
			throw new FileNotFoundException("ccggFile " + ccggFile.getAbsolutePath()
					+ " does not exist or is directory");
		}
		this.cgdFile = cgdFile;
		this.ccggFile = ccggFile;
		this.gdiFile = gdiFile;
	}
	
	public void compareToCGD()
	{
		HashMap<String, CCGGEntry> ccgg = cu.getGeneToEntry();
		HashMap<String, Integer> categoryToADARcount = new HashMap<String, Integer>();
		HashMap<String, Integer> categoryToADcount = new HashMap<String, Integer>();
		HashMap<String, Integer> categoryToARcount = new HashMap<String, Integer>();
		HashMap<String, Integer> categoryToOthercount = new HashMap<String, Integer>();
		Set<String> allCat = new HashSet<String>();
	
		for(String gene : ccgg.keySet())
		{
			if(geneToInheritance.containsKey(gene))
			{
				String cat = ccgg.get(gene).category.toString();
				String inh = geneToInheritance.get(gene);
				allCat.add(cat);
				
				if(inh.equals("AD/AR"))
				{
					if(categoryToADARcount.containsKey(cat))
					{ categoryToADARcount.put(cat, categoryToADARcount.get(cat) + 1); }
					else
					{ categoryToADARcount.put(cat, 1); }
				}
				if(inh.equals("AD"))
				{
					if(categoryToADcount.containsKey(cat))
					{ categoryToADcount.put(cat, categoryToADcount.get(cat) + 1); }
					else
					{ categoryToADcount.put(cat, 1); }
				}
				if(inh.equals("AR"))
				{
					if(categoryToARcount.containsKey(cat))
					{ categoryToARcount.put(cat, categoryToARcount.get(cat) + 1); }
					else
					{ categoryToARcount.put(cat, 1); }
				}
				if(inh.equals("Other"))
				{
					if(categoryToOthercount.containsKey(cat))
					{ categoryToOthercount.put(cat, categoryToOthercount.get(cat) + 1); }
					else
					{ categoryToOthercount.put(cat, 1); }
				}
			}
		}

		for(String cat : allCat)
		{
			System.out.print("\t" + cat);
		}
		System.out.print("\nAD/AR");
		for(String cat : allCat)
		{
			System.out.print("\t" + (categoryToADARcount.get(cat) != null ?  categoryToADARcount.get(cat) : "0"));
		}
		System.out.print("\nAD");
		for(String cat : allCat)
		{
			System.out.print("\t" + (categoryToADcount.get(cat) != null ?  categoryToADcount.get(cat) : "0"));
		}
		System.out.print("\nAR");
		for(String cat : allCat)
		{
			System.out.print("\t" + (categoryToARcount.get(cat) != null ?  categoryToARcount.get(cat) : "0"));
		}
		System.out.print("\nOther");
		for(String cat : allCat)
		{
			System.out.print("\t" + (categoryToOthercount.get(cat) != null ?  categoryToOthercount.get(cat) : "0"));
		}
		System.out.println();
		
		
		
		ArrayList<Double> popMeansCADDforADAR = new ArrayList<Double>();
		ArrayList<Double> pathoMeansCADDforADAR = new ArrayList<Double>();
		ArrayList<Double> popMeansCADDforAR = new ArrayList<Double>();
		ArrayList<Double> pathoMeansCADDforAR = new ArrayList<Double>();
		ArrayList<Double> popMeansCADDforAD = new ArrayList<Double>();
		ArrayList<Double> pathoMeansCADDforAD = new ArrayList<Double>();
		ArrayList<Double> popMeansCADDforOther = new ArrayList<Double>();
		ArrayList<Double> pathoMeansCADDforOther = new ArrayList<Double>();
		for(String gene : ccgg.keySet())
		{
			if(geneToInheritance.containsKey(gene))
			{
				Category cat = ccgg.get(gene).category;
				String inh = geneToInheritance.get(gene);
				if(cat.equals(Category.C1) || cat.equals(Category.C2) || cat.equals(Category.C3) || cat.equals(Category.C4) || cat.equals(Category.C5))
				{
					if(inh.equals("AD/AR"))
					{
						popMeansCADDforADAR.add(ccgg.get(gene).MeanPopulationCADDScore);
						pathoMeansCADDforADAR.add(ccgg.get(gene).MeanPathogenicCADDScore);
					}
					if(inh.equals("AD"))
					{
						popMeansCADDforAD.add(ccgg.get(gene).MeanPopulationCADDScore);
						pathoMeansCADDforAD.add(ccgg.get(gene).MeanPathogenicCADDScore);
					}
					if(inh.equals("AR"))
					{
						popMeansCADDforAR.add(ccgg.get(gene).MeanPopulationCADDScore);
						pathoMeansCADDforAR.add(ccgg.get(gene).MeanPathogenicCADDScore);
					}
					if(inh.equals("Other"))
					{
						popMeansCADDforOther.add(ccgg.get(gene).MeanPopulationCADDScore);
						pathoMeansCADDforOther.add(ccgg.get(gene).MeanPathogenicCADDScore);
					}
				}
			}
		}
		
		double[] popMeansCADDforADAR_ = new double[popMeansCADDforADAR.size()];
		for(int i = 0; i < popMeansCADDforADAR.size(); i++) { popMeansCADDforADAR_[i] = popMeansCADDforADAR.get(i); }
		double[] popMeansCADDforAD_ = new double[popMeansCADDforAD.size()];
		for(int i = 0; i < popMeansCADDforAD.size(); i++) { popMeansCADDforAD_[i] = popMeansCADDforAD.get(i); }
		double[] popMeansCADDforAR_ = new double[popMeansCADDforAR.size()];
		for(int i = 0; i < popMeansCADDforAR.size(); i++) { popMeansCADDforAR_[i] = popMeansCADDforAR.get(i); }
		double[] popMeansCADDforOther_ = new double[popMeansCADDforOther.size()];
		for(int i = 0; i < popMeansCADDforOther.size(); i++) { popMeansCADDforOther_[i] = popMeansCADDforOther.get(i); }
		
		double[] pathoMeansCADDforADAR_ = new double[pathoMeansCADDforADAR.size()];
		for(int i = 0; i < pathoMeansCADDforADAR.size(); i++) { pathoMeansCADDforADAR_[i] = pathoMeansCADDforADAR.get(i); }
		double[] pathoMeansCADDforAD_ = new double[pathoMeansCADDforAD.size()];
		for(int i = 0; i < pathoMeansCADDforAD.size(); i++) { pathoMeansCADDforAD_[i] = pathoMeansCADDforAD.get(i); }
		double[] pathoMeansCADDforAR_ = new double[pathoMeansCADDforAR.size()];
		for(int i = 0; i < pathoMeansCADDforAR.size(); i++) { pathoMeansCADDforAR_[i] = pathoMeansCADDforAR.get(i); }
		double[] pathoMeansCADDforOther_ = new double[pathoMeansCADDforOther.size()];
		for(int i = 0; i < pathoMeansCADDforOther.size(); i++) { pathoMeansCADDforOther_[i] = pathoMeansCADDforOther.get(i); }
		
		Mean mean = new Mean();
		double popMeanCADDforADAR = mean.evaluate(popMeansCADDforADAR_);
		double popMeanCADDforAD = mean.evaluate(popMeansCADDforAD_);
		double popMeanCADDforAR = mean.evaluate(popMeansCADDforAR_);
		double popMeanCADDforOther = mean.evaluate(popMeansCADDforOther_);
		
		double pathoMeanCADDforADAR = mean.evaluate(pathoMeansCADDforADAR_);
		double pathoMeanCADDforAD = mean.evaluate(pathoMeansCADDforAD_);
		double pathoMeanCADDforAR = mean.evaluate(pathoMeansCADDforAR_);
		double pathoMeanCADDforOther = mean.evaluate(pathoMeansCADDforOther_);
		
		System.out.println("\n\tMeanPopCADD\tMeanPathoCADD\nAD/AR\t" + popMeanCADDforADAR + "\t" + pathoMeanCADDforADAR + "\nAD\t" + popMeanCADDforAD + "\t" + pathoMeanCADDforAD+ "\nAR\t" + popMeanCADDforAR + "\t" + pathoMeanCADDforAR+ "\nOther\t" + popMeanCADDforOther + "\t" + pathoMeanCADDforOther);
		
	}
	
	public void compareToGDI()
	{
		HashMap<String, CCGGEntry> ccgg = cu.getGeneToEntry();
		
		HashMap<String, Integer> categoryToHighDmgCount = new HashMap<String, Integer>();
		HashMap<String, Integer> categoryToMediumDmgcount = new HashMap<String, Integer>();
		HashMap<String, Integer> categoryToLowDmgcount = new HashMap<String, Integer>();
		Set<String> allCat = new HashSet<String>();
		
		for(String gene : ccgg.keySet())
		{
			if(geneToDamage.containsKey(gene))
			{
				String cat = ccgg.get(gene).category.toString();
				allCat.add(cat);
				String dmg = geneToDamage.get(gene);
				if(dmg.equals("High"))
				{
					if(categoryToHighDmgCount.containsKey(cat))
					{ categoryToHighDmgCount.put(cat, categoryToHighDmgCount.get(cat) + 1); }
					else
					{ categoryToHighDmgCount.put(cat, 1); }
				}
				if(dmg.equals("Medium"))
				{
					if(categoryToMediumDmgcount.containsKey(cat))
					{ categoryToMediumDmgcount.put(cat, categoryToMediumDmgcount.get(cat) + 1); }
					else
					{ categoryToMediumDmgcount.put(cat, 1); }
				}
				if(dmg.equals("Low"))
				{
					if(categoryToLowDmgcount.containsKey(cat))
					{ categoryToLowDmgcount.put(cat, categoryToLowDmgcount.get(cat) + 1); }
					else
					{ categoryToLowDmgcount.put(cat, 1); }
				}
			}
		}
		
		System.out.println();
		for(String cat : allCat)
		{
			System.out.print("\t" + cat);
		}
		System.out.print("\nHighDmg");
		for(String cat : allCat)
		{
			System.out.print("\t" + (categoryToHighDmgCount.get(cat) != null ?  categoryToHighDmgCount.get(cat) : "0"));
		}
		System.out.print("\nMediumDmg");
		for(String cat : allCat)
		{
			System.out.print("\t" + (categoryToMediumDmgcount.get(cat) != null ?  categoryToMediumDmgcount.get(cat) : "0"));
		}
		System.out.print("\nLowDmg");
		for(String cat : allCat)
		{
			System.out.print("\t" + (categoryToLowDmgcount.get(cat) != null ?  categoryToLowDmgcount.get(cat) : "0"));
		}
		
		
		
		ArrayList<Double> popMeansCADDforHighDmg = new ArrayList<Double>();
		ArrayList<Double> pathoMeansCADDforHighDmg = new ArrayList<Double>();
		ArrayList<Double> popMeansCADDforMediumDmg = new ArrayList<Double>();
		ArrayList<Double> pathoMeansCADDforMediumDmg = new ArrayList<Double>();
		ArrayList<Double> popMeansCADDforLowDmg = new ArrayList<Double>();
		ArrayList<Double> pathoMeansCADDforLowDmg = new ArrayList<Double>();
		
		for(String gene : ccgg.keySet())
		{
			if(geneToDamage.containsKey(gene))
			{
				Category cat = ccgg.get(gene).category;
				String dmg = geneToDamage.get(gene);
				if(cat.equals(Category.C1) || cat.equals(Category.C2) || cat.equals(Category.C3) || cat.equals(Category.C4) || cat.equals(Category.C5))
				{
					if(dmg.equals("High"))
					{
						popMeansCADDforHighDmg.add(ccgg.get(gene).MeanPopulationCADDScore);
						pathoMeansCADDforHighDmg.add(ccgg.get(gene).MeanPathogenicCADDScore);
					}
					if(dmg.equals("Medium"))
					{
						popMeansCADDforMediumDmg.add(ccgg.get(gene).MeanPopulationCADDScore);
						pathoMeansCADDforMediumDmg.add(ccgg.get(gene).MeanPathogenicCADDScore);
					}
					if(dmg.equals("Low"))
					{
						popMeansCADDforLowDmg.add(ccgg.get(gene).MeanPopulationCADDScore);
						pathoMeansCADDforLowDmg.add(ccgg.get(gene).MeanPathogenicCADDScore);
					}
				}
			}
		}
		
		double[] popMeansCADDforHighDmg_ = new double[popMeansCADDforHighDmg.size()];
		for(int i = 0; i < popMeansCADDforHighDmg.size(); i++) { popMeansCADDforHighDmg_[i] = popMeansCADDforHighDmg.get(i); }
		double[] popMeansCADDforMediumDmg_ = new double[popMeansCADDforMediumDmg.size()];
		for(int i = 0; i < popMeansCADDforMediumDmg.size(); i++) { popMeansCADDforMediumDmg_[i] = popMeansCADDforMediumDmg.get(i); }
		double[] popMeansCADDforLowDmg_ = new double[popMeansCADDforLowDmg.size()];
		for(int i = 0; i < popMeansCADDforLowDmg.size(); i++) { popMeansCADDforLowDmg_[i] = popMeansCADDforLowDmg.get(i); }
		
		double[] pathoMeansCADDforHighDmg_ = new double[pathoMeansCADDforHighDmg.size()];
		for(int i = 0; i < pathoMeansCADDforHighDmg.size(); i++) { pathoMeansCADDforHighDmg_[i] = pathoMeansCADDforHighDmg.get(i); }
		double[] pathoMeansCADDforMediumDmg_ = new double[pathoMeansCADDforMediumDmg.size()];
		for(int i = 0; i < pathoMeansCADDforMediumDmg.size(); i++) { pathoMeansCADDforMediumDmg_[i] = pathoMeansCADDforMediumDmg.get(i); }
		double[] pathoMeansCADDforLowDmg_ = new double[pathoMeansCADDforLowDmg.size()];
		for(int i = 0; i < pathoMeansCADDforLowDmg.size(); i++) { pathoMeansCADDforLowDmg_[i] = pathoMeansCADDforLowDmg.get(i); }
	
		Mean mean = new Mean();
		double popMeanCADDforHighDmg = mean.evaluate(popMeansCADDforHighDmg_);
		double popMeanCADDforMediumDmg = mean.evaluate(popMeansCADDforMediumDmg_);
		double popMeanCADDforLowDmg = mean.evaluate(popMeansCADDforLowDmg_);

		double pathoMeanCADDforHighDmg = mean.evaluate(pathoMeansCADDforHighDmg_);
		double pathoMeanCADDforMediumDmg = mean.evaluate(pathoMeansCADDforMediumDmg_);
		double pathoMeanCADDforLowDmg = mean.evaluate(pathoMeansCADDforLowDmg_);
			
		System.out.println("\n\tMeanPopCADD\tMeanPathoCADD\nHighDmg\t" + popMeanCADDforHighDmg + "\t" + pathoMeanCADDforHighDmg + "\nMediumDmg\t" + popMeanCADDforMediumDmg + "\t" + pathoMeanCADDforMediumDmg + "\nLowDmg\t" + popMeanCADDforLowDmg + "\t" + pathoMeanCADDforLowDmg);
				
	}
	
	public void loadData() throws Exception
	{
		CCGGUtils cu = new CCGGUtils(ccggFile);
		
		HashMap<String, String> geneToDamage = new HashMap<String, String>();
		
		Scanner gdiScanner= new Scanner(gdiFile);

		//skip header
		gdiScanner.nextLine();
		
		//example:
		//AGPAT2	430.44157	4.33028	Medium
		while(gdiScanner.hasNextLine())
		{
			String line = gdiScanner.nextLine();
			String[] lineSplit = line.split("\t", -1);
			geneToDamage.put(lineSplit[0], lineSplit[3]);
		}
		
		
		HashMap<String, String> geneToInheritance = new HashMap<String, String>();
		HashMap<String, String> geneToOnset = new HashMap<String, String>();
		
		Scanner cgdScanner= new Scanner(cgdFile);
		
		//skip header
		cgdScanner.nextLine();
		
		//example:
		//AXIN2	Oligodontia-colorectal cancer syndrome	AD	Adult
		while(cgdScanner.hasNextLine())
		{
			String line = cgdScanner.nextLine();
			String[] lineSplit = line.split("\t", -1);
			String inheritance;
			if(lineSplit[2].contains("AD") && lineSplit[2].contains("AR"))
			{
				inheritance = "AD/AR";
			}
			else if(lineSplit[2].contains("AR"))
			{
				inheritance = "AR";
			}
			else if(lineSplit[2].contains("AD"))
			{
				inheritance = "AD";
			}
			else
			{
				inheritance = "Other";
			}
			String ageGroup;
			if(lineSplit[3].contains("Pediatric"))
			{
				ageGroup = "Pediatric";
			}
			else if(lineSplit[3].contains("Adult"))
			{
				ageGroup = "Adult";
			}
			else
			{
				ageGroup = "Other";
			}
			geneToInheritance.put(lineSplit[0], inheritance);
			geneToOnset.put(lineSplit[0], ageGroup);
			
			
		}
		
		this.cu = cu;
		this.geneToDamage = geneToDamage;
		this.geneToInheritance = geneToInheritance;
		this.geneToOnset = geneToOnset;
	}

}
