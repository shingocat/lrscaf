/*
*File: agis.ps.EdgeM5.java
*User: mqin
*Email: mqin@ymail.com
*Date: 2016年1月23日
*/
package agis.ps;

import java.io.Serializable;

import agis.ps.link.Contig;
import agis.ps.link.M5Record;
import agis.ps.util.Strand;

public class EdgeM5 implements Serializable {
	private M5Record origin;
	private M5Record terminus;
	private Strand oStrand;
	private Strand tStrand;
	private int linkNum;
	private int distMean;
	private int distSd;
	
}


