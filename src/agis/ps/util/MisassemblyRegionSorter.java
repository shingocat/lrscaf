/** 
** Usage: TODO
** Author: mqin
** Email: mqin@outlook.com
** Date: 2017年7月31日
*/
package agis.ps.util;

import java.util.Comparator;

public class MisassemblyRegionSorter implements Comparator<Object> {

	@Override
	public int compare(Object o1, Object o2) {
		MisassemblyRegion mr1 = (MisassemblyRegion) o1;
		MisassemblyRegion mr2 = (MisassemblyRegion) o2;
		if(mr1.getStart() <= mr2.getStart())
		{
			if(mr1.getEnd() <= mr2.getEnd())
				return 1;
			else 
				return 0;
		} else
		{
			return 0;
		}
	}

}
