package org.molgenis.calibratecadd.structs;

import java.util.ArrayList;
import java.util.HashMap;

public class ChrPosRefAltUniqueVariants
{

	//VCF lines
	private HashMap<String,String> lines;
	
	//duplicates..
	private ArrayList<String> duplicateLines;
	
	public ChrPosRefAltUniqueVariants()
	{
		this.lines = new HashMap<String, String>();
		this.duplicateLines = new ArrayList<String>();
	}
	
	public boolean add(String line)
	{
		String[] split = line.split("\t", -1);
		String key = split[0] + "_" + split[1] + "_" + split[3] + "_" + split[4];
		if(!lines.containsKey(key))
		{
			lines.put(key, line);
			return true;
		}
		else
		{
			duplicateLines.add(line);
			return false;
		}
	}
	
	public HashMap<String,String> getLines()
	{
		return lines;
	}
	
	public ArrayList<String> getDuplicateLines()
	{
		return duplicateLines;
	}
}
