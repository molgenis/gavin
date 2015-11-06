package org.molgenis.calibratecadd.support;

import java.util.HashMap;
import java.util.List;

import org.molgenis.data.Entity;

public class EntityPlus
{

	private Entity e;
	private HashMap<String, Object> keyVal;

	public EntityPlus(Entity e)
	{
		this.e = e;
		this.keyVal = new HashMap<String, Object>();
	}

	public HashMap<String, Object> getKeyVal()
	{
		return keyVal;
	}

	public void setKeyVal(HashMap<String, Object> keyVal)
	{
		this.keyVal = keyVal;
	}

	public Entity getE()
	{
		return e;
	}

	@Override
	public String toString()
	{
		return "EntityPlus [e=" + e + ", keyVal=" + keyVal + "]";
	}
	

	public static boolean contains(List<EntityPlus> list, Entity ep)
	{
		for(EntityPlus e : list)
		{
			if(e.getE().equals(ep))
			{
				return true;
			}
		}
		return false;
	}

}
