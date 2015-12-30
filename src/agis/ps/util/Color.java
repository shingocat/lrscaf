/*
*File: agis.ps.util.Color.java
*User: mqin
*Email: mqin@ymail.com
*Date: 2015Äê12ÔÂ29ÈÕ
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


