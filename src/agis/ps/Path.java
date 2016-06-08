/*
*File: agis.ps.Path.java
*User: mqin
*Email: mqin@ymail.com
*Date: 2015年12月28日
*/
package agis.ps;

import java.util.LinkedList;

import agis.ps.link.Edge;
import agis.ps.seqs.Contig;
import agis.ps.util.Strand;

public class Path {
	// the vertices and verStatus should be synchronized, i.e. add or remove.
	private LinkedList<Contig> vertices; // store contig;
	private LinkedList<Strand> verStatus; // store strand over each contig;
	
	public Path()
	{
		if(vertices == null)
			vertices = new LinkedList<Contig>();
		vertices.clear();
		if(verStatus == null)
			verStatus = new LinkedList<Strand>();
		verStatus.clear();
	}
	
	//push one vertex Strand status in the verStatus linked list;
	public void pushStrand(Strand s)
	{
		verStatus.addLast(s);
	}
	
	// unshift one vertex strand status in the verStatus linked list;
	public void unshiftStrand(Strand s)
	{
		verStatus.addFirst(s);
	}
	
	public Strand getStrandByIndex(int index)
	{
		if(verStatus == null || verStatus.isEmpty())
			return null;
		return verStatus.get(index);
	}
	
	// push one vertex to the end of this linked list;
	public void push(Contig v)
	{
		vertices.addLast(v);
	}
	
	// pop one vertex from the end of this linked list;
	// if the linked list is empty, do nothing;
	public Contig pop()
	{
		return vertices.pollLast();
	}
	
	// shift one vertex from the first of this linked list;
	public Contig shift()
	{
		return vertices.pollFirst();
	}
	
	// unshift one vertex to the first of this linked list;
	public void unshift(Contig v)
	{
		vertices.addFirst(v);
	}
	
	// get the first element 
	public Contig getFirstElement()
	{
		return vertices.getFirst();
	}
	
	// get the last element
	public Contig getLastElement()
	{
		return vertices.getLast();
	}
	
	// equal the first element
	public boolean isEqualFirstElement(Contig c)
	{
		Contig f = this.getFirstElement();
		if(f.getID().equals(c.getID()))
			return true;
		return false;
	}
	//equal the last element
	public boolean isEqualLastElement(Contig c)
	{
		Contig f = this.getLastElement();
		if(f.getID().equals(c.getID()))
			return true;
		return false;
	}
	// path cotain reverse edge, i.e. A->B->C in the path, and now edge C->B is the reverse edge in the path
	public boolean isExistReverseEdge(Edge e)
	{
		boolean value = false;
		String oId = e.getOrigin().getID();
		String tId = e.getTerminus().getID();
		for(int i = getSize() - 1; i > 0; i--)
		{
			if(this.vertices.get(i).getID().equals(oId) && this.vertices.get(i - 1).getID().equals(tId))
			{
				value = true;
				break;
			}
		}
		return value;
	}
	
	// path cotain edge, i.e. A->B->C int the path, and now edge B->C is the edge;
	public boolean isExistEdge(Edge e)
	{
		boolean value = false;
		String oId = e.getOrigin().getID();
		String tId = e.getTerminus().getID();
		for(int i = 0; i < getSize() - 1; i++)
		{
			if(this.vertices.get(i).getID().equals(oId) && this.vertices.get(i + 1).getID().equals(tId))
			{
				value = true;
				break;
			}
		}
		return value;
	}
	
	// get path element by index
	public Contig getElement(int index)
	{
		if(index > getSize() - 1)
		{
			throw new IllegalArgumentException("The index out of box in path, larger than size!");
		}
		return this.vertices.get(index);
	}
	
	public int getSize()
	{
		return this.vertices.size();
	}
	
	public boolean isEmpty()
	{
		if(getSize() == 0)
			return true;
		return false;
	}
	
	//return -1 means that do not contain this contig in path, else return the index of this contig in the path
	public int containVertex(Contig origin) {
		boolean isExist = getVertices().contains(origin);
		if(isExist)
		{
			int index = 0;
			for(Contig c : getVertices())
			{
				if(c.getID().equals(origin.getID()))
				{
					return index;
				}
				index++;
			}
		}
		return -1;
	}
	
	public LinkedList<Contig> getVertices()
	{
		return this.vertices;
	}
	
	@Override
	public String toString()
	{
		StringBuilder sb = new StringBuilder();
		for(int i = 0 ; i < getSize(); i++)
		{
			sb.append(getElement(i).getID());
			if(i != getSize() - 1)
				sb.append("->");
		}
		return sb.toString();
	}
}


