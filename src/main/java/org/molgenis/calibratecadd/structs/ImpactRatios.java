package org.molgenis.calibratecadd.structs;

public class ImpactRatios
{

	public double high;
	public double moderate;
	public double low;
	public double modifier;
	
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
		//if we have 33+33+33+0 or so, add 1 to 'moderate'
		if(Math.round(total) != 100)
		{
			throw new Exception("Impact rations should add up to ~100, but instead: " + high + "+" + moderate + "+" + low + "+"+modifier +"=" + (high+moderate+low+modifier));
		}
	}

	@Override
	public String toString()
	{
		return "ImpactRatios [high=" + high + ", moderate=" + moderate + ", low=" + low + ", modifier=" + modifier
				+ "]";
	}
	
	
	
	
}
