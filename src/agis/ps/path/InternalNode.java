/*
*File: agis.ps.path.InternalPath2.java
*User: mqin
*Email: mqin@ymail.com
*Date: 2017年1月17日
*/
package agis.ps.path;

import agis.ps.seqs.Contig;

public class InternalNode {
	private InternalNode parent;
	private Contig child;
	private boolean isLeaf;
	public InternalNode getParent() {
		return parent;
	}
	public void setParent(InternalNode parent) {
		this.parent = parent;
	}
	public Contig getChild() {
		return child;
	}
	public void setChild(Contig child) {
		this.child = child;
	}
	public boolean isLeaf() {
		return isLeaf;
	}
	public void setLeaf(boolean isLeaf) {
		this.isLeaf = isLeaf;
	}
	@Override
	public String toString() {
		return "InternalPath2 [parent=" + parent + ", chile=" + child + ", isLeaf=" + isLeaf + "]";
	}	
}


