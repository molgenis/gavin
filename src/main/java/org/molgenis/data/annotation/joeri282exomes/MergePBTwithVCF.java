package org.molgenis.data.annotation.joeri282exomes;

import java.io.File;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;
import java.util.Scanner;
import java.util.stream.Collectors;

import org.molgenis.data.Entity;
import org.molgenis.data.annotation.cmd.CommandLineAnnotatorConfig;
import org.molgenis.data.vcf.VcfRepository;
import org.molgenis.data.vcf.utils.VcfUtils;
import org.springframework.context.annotation.AnnotationConfigApplicationContext;
import org.springframework.stereotype.Component;

@Component
public class MergePBTwithVCF
{

	HashMap<String, ArrayList<String>> patientgroupToGenes = new HashMap<String, ArrayList<String>>();

	private ArrayList<String> pbtChrPos = new ArrayList<String>();

	public void go(File vcfFile, File pbtFile, File output) throws Exception
	{

		PrintWriter pw = new PrintWriter(output, "UTF-8");

		Scanner s = new Scanner(pbtFile);
		String line = null;
		while (s.hasNextLine())
		{
			line = s.nextLine();
			String[] lineSplit = line.split("\t", -1);
			pbtChrPos.add(lineSplit[0] + "_" + lineSplit[1]);
		}
		s.close();
		
		
		
		Map<String, Long> counts = pbtChrPos.stream().collect(Collectors.groupingBy(e -> e, Collectors.counting()));
		
		

		Iterator<Entity> vcf = new VcfRepository(vcfFile, "vcf").iterator();

		while (vcf.hasNext())
		{
			Entity record = vcf.next();

			String chr = record.getString("#CHROM");
			String pos = record.getString("POS");

			
			if (counts.containsKey(chr + "_" + pos))
			{
				for(int i = 0; i < counts.get(chr + "_" + pos); i++)
				{
					pw.println(VcfUtils.convertToVCF(record));
					
				}
				pw.flush();
			}

		}
		pw.flush();
		pw.close();
	}

	public static void main(String[] args) throws Exception
	{
		// configureLogging();

		// See http://stackoverflow.com/questions/4787719/spring-console-application-configured-using-annotations
		AnnotationConfigApplicationContext ctx = new AnnotationConfigApplicationContext("org.molgenis.data.annotation");
		ctx.register(CommandLineAnnotatorConfig.class);
		MergePBTwithVCF main = ctx.getBean(MergePBTwithVCF.class);
		main.run(args);
		ctx.close();
	}

	public void run(String[] args) throws Exception
	{
		if (!(args.length == 3))
		{
			throw new Exception("Must supply 3 arguments");
		}

		File vcfFile = new File(args[0]);
		if (!vcfFile.isFile())
		{
			throw new Exception("Input VCF file does not exist or directory: " + vcfFile.getAbsolutePath());
		}

		File pbtFile = new File(args[1]);
		if (!pbtFile.isFile())
		{
			throw new Exception("PBT file does not exist or directory: " + pbtFile.getAbsolutePath());
		}

		File output = new File(args[2]);
		if (output.exists())
		{
			System.out.println("WARNING: output file already exists, overwriting " + output);
		}

		MergePBTwithVCF cf = new MergePBTwithVCF();
		cf.go(vcfFile, pbtFile, output);

	}

}
