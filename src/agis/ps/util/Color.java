/*
*File: agis.ps.util.Color.java
*User: mqin
*Email: mqin@ymail.com
*Date: 2015年12月28日
*/
package agis.ps.util;

public enum Color {
	RED("red"), BLUE("blue"), GREEN("green"), BROWN("brown");
	
	private String type;
	
	private Color(String type)
	{
		this.type = type;
	}
	
	public String toString()
	{
		return type;
	}
}


