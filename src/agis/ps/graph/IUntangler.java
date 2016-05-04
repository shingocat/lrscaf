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

import agis.ps.Edge;
import agis.ps.link.Contig;

public interface IUntangler {
	final static Logger logger = LoggerFactory.getLogger(IUntangler.class);
	
	public void transitiveReducting();
	
	public void linearMergin();
}


