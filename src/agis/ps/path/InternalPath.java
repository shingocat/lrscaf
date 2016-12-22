/*
*File: agis.ps.path.InternalPath.java
*User: mqin
*Email: mqin@ymail.com
*Date: 2016年12月22日
*/
package agis.ps.path;

import java.util.LinkedList;

import agis.ps.seqs.Contig;

public class InternalPath {
	private LinkedList<Contig> path = new LinkedList<Contig>(); // the Linked contig list for this internal path
	private int score = 0; // this path score
	
	public void addFirst(Contig fst)
	{
		path.addFirst(fst);
	}
	
	public void addLast(Contig lst)
	{
		path.addLast(lst);
	}
	
	public void addScore(int score)
	{
		this.score += score;
	}
	
	public void subtractScore(int score)
	{
		this.score -= score;
	}
	
	public int getScore()
	{
		return this.score;
	}
	
	public boolean isContain(Contig c)
	{
		if(path.contains(c))
			return true;
		else 
			return false;
	}
	
	public int getIndex(Contig c)
	{
		return path.indexOf(c);
	}
	
	public boolean isEmpty()
	{
		if(path == null || path.isEmpty())
			return true;
		else 
			return false;
	}
	
	public LinkedList<Contig> getPath()
	{
		return path;
	}

	@Override
	public String toString() {
		return "InternalPath [path=" + path + ", score=" + score + "]";
	}
}


