/** 
** Usage: TODO
** Author: mqin
** Email: mqin@outlook.com
** Date: 2019年1月23日
*/
package agis.ps.link;

import java.util.Comparator;

public class ByLocOrderComparator implements Comparator<Object> {
	
	@Override
	public int compare(Object cFirst, Object cSecond) {
		MRecord m1 = (MRecord) cFirst;
		MRecord m2 = (MRecord) cSecond;
		int diff = Integer.valueOf(m1.getqStart()) - Integer.valueOf(m2.getqStart());
		if (diff > 0)
			return 1;
		if (diff < 0)
			return -1;
		else
			return 0;
	}
	
}
