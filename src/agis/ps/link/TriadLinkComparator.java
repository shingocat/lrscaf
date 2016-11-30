/*
*File: agis.ps.link.TriadLinkComparator.java
*User: mqin
*Email: mqin@ymail.com
*Date: 2016年5月14日
*/
package agis.ps.link;

import java.util.Comparator;

public class TriadLinkComparator implements Comparator<TriadLink> {

	@Override
	public int compare(TriadLink o1, TriadLink o2) {
		// TODO Auto-generated method stub
		if(o1.getSupLinks() > o2.getSupLinks())
			return -1;
		else if(o1.getSupLinks() < o2.getSupLinks())
			return 1;
		else
			return 0;
	}

}


