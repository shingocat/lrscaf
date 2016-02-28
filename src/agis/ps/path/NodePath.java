/*
*File: agis.ps.path.NodePath.java
*User: mqin
*Email: mqin@ymail.com
*Date: 2016年2月26日
*/
package agis.ps.path;

import java.io.Serializable;
import java.util.LinkedList;

import agis.ps.exception.EmptyPathException;
import agis.ps.link.Contig;


public class NodePath implements Serializable{
		
	private LinkedList<Node> path = new LinkedList<Node>();
	
	public void push(Node node)
	{
		path.addLast(node);
	}
	
	public Node pop()
	{
		return path.pollLast();
	}
	
	public Node shift()
	{
		return path.pollFirst();
	}
	
	public void unshift(Node node)
	{
		path.addFirst(node);
	}
	
	public int getPathSize()
	{
		return path.size();
	}
	
	// checking whether the contig in the path by their index;
	// index == 0, denoted as checking all the element of path,
	// index == positive number, denoted as checking the forward specified index, 
	// such as index == 1, means checking the first element,
	// index == negative number, denoted as checking the reverse specified index,
	// such as index == -1, means, checking the last element,
	public boolean isNextExist(Contig next, int index) throws EmptyPathException {
		// TODO Auto-generated method stub
		if(path == null)
			path = new LinkedList<Node>();
		if(path.isEmpty())
			return false;
		boolean isExist = false;
		if(index == 0)
		{
			for(int i =0 ; i < path.size(); i++)
			{
				Node node = path.get(i);
				if(node.getCnt().equals(next))
				{	
					isExist = true;
					break;
				}
			}
		} else if (index < 0)
		{
			int pathSize = path.size();
			if(pathSize + index < 0)
			return false;
//				throw new IllegalArgumentException("The index is out of negative boundary!");
			Node node = path.get(pathSize + index);
			if(node.getCnt().equals(next))
				isExist = true;
		} else if (index > 0)
		{
			Node node = path.get(index - 1);
			if(node.getCnt().equals(next))
				isExist = true;
		}
		return isExist;
	}
	
	public Node getElement(int index)
	{
		return path.get(index);
	}
	
	@Override
	public String toString()
	{
		StringBuilder sb = new StringBuilder();
		for(int i = 0 ; i < getPathSize(); i++)
		{
			sb.append(getElement(i).getCnt().getID());
			if(i != getPathSize() - 1)
				sb.append("->");
		}
		return sb.toString();
	}
}


