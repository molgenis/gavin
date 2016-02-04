package org.molgenis.data.annotation.joeri282exomes.struct;

/**
 * 
 * We store counts for acting genotype counts for variants in 3 types : dominant, additive and recessive.
 * In addition, we also store the 'non acting' genotypes under each of those models.
 * 
 * Counts for acting:
 * 		0/0	0/1	1/1
 * Dom.	0	1	1
 * Add.	0	1	2
 * Rec.	0	0	1
 * 
 * Counts for non-acting (=inverse):
 * 		0/0	0/1	1/1
 * Dom.	1	0	0
 * Add.	2	1	0
 * Rec.	1	1	0
 * 
 * Additive is on a different scale than Dom/Rec to account for the partial heterozygous effect.
 * 
 * Order: acting dominant, acting additive, acting recessive | non-acting dominant, non-acting additive, non-acting recessive
 * 
 * toString examples:
 * 3,3,0|100,500,400
 * 5,8,3|200,750,550
 * 
 * @author jvelde
 *
 */
public class AlleleCounts
{

	public static String printNull = "0,0,0|0,0,0";
	
	private int actingDominant;
	private int actingAdditive;
	private int actingRecessive;
	private int nonActingDominant;
	private int nonActingAdditive;
	private int nonActingRecessive;
	
	public AlleleCounts()
	{
		super();
		this.actingDominant = 0;
		this.actingAdditive = 0;
		this.actingRecessive = 0;
		this.nonActingDominant = 0;
		this.nonActingAdditive = 0;
		this.nonActingRecessive = 0;
	}
	
	public AlleleCounts(String acString) throws Exception
	{
		String[] actingVsNonActing = acString.split("\\|", -1);
		if(actingVsNonActing.length != 2)
		{
			throw new Exception("acting vs non acting split != 2");
		}
		
		String[] acting = actingVsNonActing[0].split(",", -1);
		String[] nonActing = actingVsNonActing[1].split(",", -1);
		
		if(acting.length != 3 || nonActing.length != 3)
		{
			throw new Exception("acting OR non acting split != 3");
		}
		
		this.actingDominant = Integer.parseInt(acting[0]);
		this.actingAdditive = Integer.parseInt(acting[1]);
		this.actingRecessive = Integer.parseInt(acting[2]);
		this.nonActingDominant = Integer.parseInt(nonActing[0]);
		this.nonActingAdditive = Integer.parseInt(nonActing[1]);
		this.nonActingRecessive = Integer.parseInt(nonActing[2]);
	}

	public void add(AlleleCounts ac)
	{
		this.actingDominant += ac.getActingDominant();
		this.actingAdditive += ac.getActingAdditive();
		this.actingRecessive += ac.getActingRecessive();
		this.nonActingDominant += ac.getNonActingDominant();
		this.nonActingAdditive += ac.getNonActingAdditive();
		this.nonActingRecessive += ac.getNonActingRecessive();
	}
	
	

	public int getActingDominant()
	{
		return actingDominant;
	}



	public void setActingDominant(int actingDominant)
	{
		this.actingDominant = actingDominant;
	}



	public int getActingAdditive()
	{
		return actingAdditive;
	}



	public void setActingAdditive(int actingAdditive)
	{
		this.actingAdditive = actingAdditive;
	}



	public int getActingRecessive()
	{
		return actingRecessive;
	}



	public void setActingRecessive(int actingRecessive)
	{
		this.actingRecessive = actingRecessive;
	}



	public int getNonActingDominant()
	{
		return nonActingDominant;
	}



	public void setNonActingDominant(int nonActingDominant)
	{
		this.nonActingDominant = nonActingDominant;
	}



	public int getNonActingAdditive()
	{
		return nonActingAdditive;
	}



	public void setNonActingAdditive(int nonActingAdditive)
	{
		this.nonActingAdditive = nonActingAdditive;
	}



	public int getNonActingRecessive()
	{
		return nonActingRecessive;
	}



	public void setNonActingRecessive(int nonActingRecessive)
	{
		this.nonActingRecessive = nonActingRecessive;
	}



	@Override
	public String toString()
	{
		return actingDominant+","+actingAdditive+","+actingRecessive+"|"+nonActingDominant+","+nonActingAdditive+","+nonActingRecessive;
	}
	
	

}
