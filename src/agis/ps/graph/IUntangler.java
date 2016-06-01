/*
*File: agis.ps.util.Untangler.java
*User: mqin
*Email: mqin@ymail.com
*Date: 2016年4月27日
*/
package agis.ps.graph;

import java.util.List;
import java.util.Vector;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import agis.ps.link.Contig;
import agis.ps.link.Edge;

public interface IUntangler {
	final static Logger logger = LoggerFactory.getLogger(IUntangler.class);
	
	public void transitiveReducting();
	
	public void linearMergin();
	// delete the error prone edge by ratio, it mean if one node have more than 2 adjacent vertex, if the
	// edge A->B(10), C->B(20) and D->B(2), the edge D->B(2) might be error, so delete it;
	public void delErrorProneEdge(double ratio);
}


