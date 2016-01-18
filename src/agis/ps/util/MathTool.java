/*
*File: agis.ps.util.MathTool.java
*User: mqin
*Email: mqin@ymail.com
*Date: 2016年1月13日
*/
package agis.ps.util;

import java.util.Arrays;
import java.util.List;

public class MathTool {
	
	public static Integer mean(List<Integer> nums)
	{
		return MathTool.mean(nums.toArray(new Integer[nums.size()]));
	}
	
	public static Integer mean(Integer [] nums)
	{
		int sum = 0;
		int size = nums.length;
		for(int i = 0; i < size; i++)
		{
			sum += nums[i];
		}
		return Math.round(sum / size);
	}
	
	public static Integer sd(List<Integer> nums)
	{
		return MathTool.sd(nums.toArray(new Integer[nums.size()]));
	}
	
	public static Integer sd(Integer [] nums)
	{
		int mean = MathTool.mean(nums);
		int size = nums.length;
		double diff = 0.0f;
		for(int i = 0; i < size; i++)
		{
			diff += Math.pow(nums[i] - mean, 2);		
		}
		
		return (int) Math.round(Math.sqrt(diff/(size - 1)));
	}

	public static int max(int typeA, int typeB, int typeC, int typeD) {
		// TODO Auto-generated method stub
		int [] arr = {typeA, typeB, typeC, typeD};
		Arrays.sort(arr);
		return arr[3];
	}

	public static int min(Integer[] values) {
		// TODO Auto-generated method stub
		
		return 0;
	}
	
	public static int min(List<Integer> values)
	{
		return 0;
	}

}


