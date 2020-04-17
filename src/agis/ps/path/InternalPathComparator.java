/*
*File: agis.ps.path.InternalPathComparator.java
*User: mqin
*Email: mqin@ymail.com
*Date: 2016年12月22日
*/
package agis.ps.path;

import java.util.Comparator;

import agis.ps.path.InternalPath;

public class InternalPathComparator implements Comparator<InternalPath> {
	@Override
	public int compare(InternalPath o1, InternalPath o2) {
		// TODO Auto-generated method stub
		if(o1.getScore() > o2.getScore())
			return -1;
		else if(o1.getScore() < o2.getScore())
			return 1;
		else
			return 0;
	}
}


