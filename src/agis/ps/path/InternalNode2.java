/*
*File: agis.ps.path.InternalNode.java
*User: mqin
*Email: mqin@ymail.com
*Date: 2016年12月22日
*/
package agis.ps.path;

import agis.ps.seqs.Contig;

public class InternalNode2 {
	private Contig grandfather = null;
	private Contig parent = null;
	private Contig children = null;
	private boolean isLeaf = false;
	
	public void setGrandfather(Contig grandfather)
	{
		this.grandfather = grandfather;
	}
	
	public Contig getGrandfather()
	{
		return this.grandfather;
	}
	
	public void setChildren(Contig children)
	{
		this.children = children;
	}
	
	public Contig getChildren()
	{
		return this.children;
	}
	
	public void setParent(Contig parent)
	{
		this.parent = parent;
	}
	
	public Contig getParent()
	{
		return parent;
	}
	
	public boolean isLeaf()
	{
		return isLeaf;
	}
	
	public void setLeaf(boolean isLeaf)
	{
		this.isLeaf = isLeaf;
	}
	
	@Override
	public String toString() {
		return "[Grandfather=" + grandfather + ", parent=" + parent + ", children=" + children + ", isLeaf=" + isLeaf + "]";
	}
}


