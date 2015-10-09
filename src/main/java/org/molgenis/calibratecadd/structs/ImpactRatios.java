package org.molgenis.calibratecadd.structs;

public class ImpactRatios
{

	public int high;
	public int moderate;
	public int low;
	public int modifier;
	
	public ImpactRatios(int high, int moderate, int low, int modifier) throws Exception
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
		int total = high+moderate+low+modifier;
		//if we have 33+33+33+0 or so, add 1 to 'moderate'
		if(total == 99)
		{
			moderate++;
		}
		else if(total == 101)
		{
			if(moderate > 0)
			{
				moderate--;
			}
			else
			{
				throw new Exception("total is 101 and could not fix by subtracting from moderate");
			}
		}
		else if(total != 100)
		{
			throw new Exception("Impact rations should add up to 100, but instead: " + high + "+" + moderate + "+" + low + "+"+modifier +"=" + (high+moderate+low+modifier));
		}
	}

	@Override
	public String toString()
	{
		return "ImpactRatios [high=" + high + ", moderate=" + moderate + ", low=" + low + ", modifier=" + modifier
				+ "]";
	}
	
	
	
	
}
