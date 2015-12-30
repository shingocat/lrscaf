/*
*File: agis.ps.SimplePath.java
*User: mqin
*Email: mqin@ymail.com
*Date: 2015Äê12ÔÂ29ÈÕ
*/
package agis.ps;

import agis.ps.util.Color;

public class SimplePath extends Path{
	private String start;
	private String end;
	private Color color;
	private String label;
	
	public SimplePath()
	{
		
	}
	
	public SimplePath(String start, String end, Color color, String label)
	{
		this.start = start;
		this.end = end;
		this.color = color;
		this.label = label;
	}

	public String getStart() {
		return start;
	}

	public void setStart(String start) {
		this.start = start;
	}

	public String getEnd() {
		return end;
	}

	public void setEnd(String end) {
		this.end = end;
	}

	public Color getColor() {
		return color;
	}

	public void setColor(Color color) {
		this.color = color;
	}

	public String getLabel() {
		return label;
	}

	public void setLabel(String label) {
		this.label = label;
	}
	
	

}


