package org.molgenis.calibratecadd.structs;

import java.text.DecimalFormat;
import java.text.NumberFormat;

public class ImpactRatios
{

	public double high;
	public double moderate;
	public double low;
	public double modifier;
	
	NumberFormat f = new DecimalFormat("#0.00");  
	
	public ImpactRatios(double high, double moderate, double low, double modifier) throws Exception
	{
		super();
		this.high = high;
		this.moderate = moderate;
		this.low = low;
		this.modifier = modifier;
		check();
	}
	
	public void check() throws Exception
	{
		double total = high+moderate+low+modifier;
		if(Math.round(total) != 100)
		{
			throw new Exception("Impact rations should add up to ~100, but instead: " + high + "+" + moderate + "+" + low + "+"+modifier +"=" + (high+moderate+low+modifier));
		}
	}

	@Override
	public String toString()
	{
		return "high: " + f.format(high) + "%, moderate: " + f.format(moderate) + "%, low: " + f.format(low) + "%, modifier: " + f.format(modifier)+ "%";
	}
	
	
	
	
}
