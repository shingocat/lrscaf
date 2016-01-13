/*
*File: agis.ps.util.ContigPairType.java
*User: mqin
*Email: mqin@ymail.com
*Date: 2016年1月13日
*/
package agis.ps.util;

public enum ContigPairType {
	PP("+","+"), PM("+","-"),MM("-","-"),MP("-","+");
	
	String origin;
	String terminus;
	
	ContigPairType(String origin, String terminus)
	{
		this.origin  = origin;
		this.terminus = terminus;
	}
	
	public String getOrigin()
	{
		return origin;
	}
	
	public String getTerminus()
	{
		return terminus;
	}
}	


