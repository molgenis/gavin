package org.molgenis.data.annotation.joeri282exomes;

import java.io.File;
import java.io.IOException;
import java.util.Iterator;
import java.util.Map;

import org.molgenis.data.Entity;
import org.molgenis.data.annotation.RepositoryAnnotator;
import org.molgenis.data.annotation.cmd.CommandLineAnnotatorConfig;
import org.molgenis.data.vcf.VcfRepository;
import org.molgenis.util.ApplicationContextProvider;
import org.springframework.context.ApplicationContext;
import org.springframework.context.annotation.AnnotationConfigApplicationContext;
import org.springframework.stereotype.Component;

@Component
public class ClinicalFilters
{
	

	
	public void bla(File vcfFile, File patientGroups, File exacFile) throws IOException
	{
		VcfRepository vcf = new VcfRepository(vcfFile, "vcf");
		ApplicationContext applicationContext = ApplicationContextProvider.getApplicationContext();
		Map<String, RepositoryAnnotator> annotators = applicationContext.getBeansOfType(RepositoryAnnotator.class);
		RepositoryAnnotator exacAnnotator = annotators.get("exac");
		
		exacAnnotator.getCmdLineAnnotatorSettingsConfigurer().addSettings(exacFile.getAbsolutePath());
		//annotate(exacAnnotator, inputVcfFile, outputVCFFile, attrNames);
	

		
		Iterator<Entity> it = exacAnnotator.annotate(vcf);
		
		while(it.hasNext())
		{
			Entity record = it.next();
			System.out.println(record.toString());
		}
	}
	
	
	public static void main(String[] args) throws Exception
	{
	//	configureLogging();

		// See http://stackoverflow.com/questions/4787719/spring-console-application-configured-using-annotations
		AnnotationConfigApplicationContext ctx = new AnnotationConfigApplicationContext("org.molgenis.data.annotation");
		ctx.register(CommandLineAnnotatorConfig.class);
		ClinicalFilters main = ctx.getBean(ClinicalFilters.class);
		main.run(args);
		ctx.close();
	}
	

	
	public void run(String[] args) throws Exception
	{
		if (!(args.length == 3))
		{
			throw new Exception(
					"Must supply 3 arguments");
		}
		
		File vcfFile = new File(args[0]);
		if (!vcfFile.isFile())
		{
			throw new Exception("Input VCF file does not exist or directory: " + vcfFile.getAbsolutePath());
		}
		
		File patientGroups = new File(args[1]);
		if (!patientGroups.isFile())
		{
			throw new Exception("Patient groups file does not exist or directory: " + patientGroups.getAbsolutePath());
		}
		
		File exacFile = new File(args[2]);
		if (!exacFile.isFile())
		{
			throw new Exception("Exac file does not exists at " + exacFile);
		}

		ClinicalFilters cf = new ClinicalFilters();
		cf.bla(vcfFile, patientGroups, exacFile);
		
		
	}
	
}
