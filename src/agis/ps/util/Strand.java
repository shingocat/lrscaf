/*
*File: agis.ps.util.Strand.java
*User: mqin
*Email: mqin@ymail.com
*Date: 2015年12月28日
*/
package agis.ps.util;

public enum Strand {
	FORWARD("+"), REVERSE("-");
	
	private String type;
	
	Strand(String type)
	{
		this.type = type;
	}
	
	public String getType()
	{
		return type;
	}
	
	@Override
	public String toString()
	{
		return type;
	}
	
}


