package org.molgenis.calibratecadd.misc;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.Writer;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Scanner;

import freemarker.template.Configuration;
import freemarker.template.Template;
import freemarker.template.TemplateException;

public class CreateWebsite
{
	/**
	 * Arguments: [0] the output of Step 7, e.g. E:\Data\clinvarcadd\clinvar.patho.fix.snpeff.exac.genesumm.tsv [1] HTML
	 * file to write to, e.g. E:\Data\clinvarcadd\index.html
	 * 
	 * @param args
	 * @throws IOException
	 * @throws TemplateException
	 */
	public static void main(String[] args) throws IOException, TemplateException
	{

		Scanner s = new Scanner(new File(args[0]));
		s.nextLine(); //skip header
		ArrayList<String[]> genes = new ArrayList<String[]>();
		String line = null;
		while(s.hasNextLine())
		{
			line = s.nextLine();
			genes.add(line.split("\t", -1));
		}
		s.close();
		
		// Freemarker configuration object
		Configuration cfg = new Configuration();

		// Load template from source folder
		Template template = cfg.getTemplate("src/main/java/org/molgenis/calibratecadd/misc/website.ftl");

		// Build the data-model
		Map<String, Object> data = new HashMap<String, Object>();
		data.put("genes", genes);

		// File output
		Writer file = new FileWriter(new File(args[1]));
		template.process(data, file);
		file.flush();
		file.close();

	}
}