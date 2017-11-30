/*
*File: agis.ps.util.MathTool.java
*User: mqin
*Email: mqin@ymail.com
*Date: 2016年1月13日
*/
package agis.ps.util;

import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
//import java.util.Vector;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

public class MathTool {
	private static Logger logger = LoggerFactory.getLogger(MathTool.class);

	public static Long sum(List<Integer> nums) {
		long sum = 0;
		for (int i = 0; i < nums.size(); i++) {
			sum += nums.get(i);
		}
		return sum;
	}
	
	public static int sum(int [] nums)
	{
		int sum = 0;
		for (int i = 0; i < nums.length; i++) {
			sum += nums[i];
		}
		return sum;
	}
	
	public static Map<String, Double> summary(List<Integer> nums)
	{
		return MathTool.summary(nums.toArray(new Integer[nums.size()]));
	}
	
	public static Map<String, Double> summary(Integer [] nums)
	{
		Map<String, Double> values = new HashMap<String, Double>();
		int len = nums.length;
		double min = 0;
		double firstQ = 0;
		double median = 0;
		double thirdQ = 0;
		double max = 0;
		Arrays.sort(nums);
		if(len != 0)
		{
			if(len % 2 == 1)
			{ // odd number
				firstQ = median(Arrays.copyOfRange(nums, 0, len / 2 + 1), true);
				thirdQ = median(Arrays.copyOfRange(nums, len/2, len), true);
			} else
			{ // even number
				firstQ = median(Arrays.copyOfRange(nums, 0, len /2), true);
				thirdQ = median(Arrays.copyOfRange(nums, len /2, len), true);
			}
			median = median(nums, true);
			min = nums[0];
			max = nums[len - 1];
		}
		values.put("MIN", min);
		values.put("FIRSTQ", firstQ);
		values.put("MEDIAN", median);
		values.put("THIRDQ", thirdQ);
		values.put("MAX", max);
		return values;
	}	
	
	public static double median(List<Integer> nums, boolean isSorted)
	{
		return MathTool.median(nums.toArray(new Integer[nums.size()]), isSorted);
	}
	
	public static double median(Integer [] nums, boolean isSorted)
	{
		double median = 0;
		int len = nums.length;
		if(len == 0)
			return median;
		if(!isSorted)
			Arrays.sort(nums);
		if((len % 2) == 1)
		{ // odd number
			int index = (len + 1)/ 2;
			index = index - 1; // zero-based index;
			median = nums[index];
		} else
		{ // even nubmer
			int n = len/2;
			int f = n + 1;
			// zero-based index
			n = n - 1;
			f = f - 1;
			int sum = nums[n] + nums[f];
			median = (double)sum / 2;
		}
		return median;
	}
	
	public static Integer mean(List<Integer> nums) {
		return MathTool.mean(nums.toArray(new Integer[nums.size()]));
	}

	public static Integer mean(Integer [] nums) {
		int mean = 0;
		try {
			long sum = 0;
			int size = nums.length;
			if (size == 0)
				return 0;
//				throw new Exception("The arrays could not be empty when computed mean!");
			for (int i = 0; i < size; i++) {
				sum += nums[i];
			}
			mean = Math.round(sum / size);
		} catch (Exception e) {
			logger.error(MathTool.class.getName() + "\t" + e.getMessage());
		}
		return mean;
	}
	
	public static int mean(int [] nums)
	{
		int mean = 0;
		try {
			long sum = 0;
			int size = nums.length;
			if (size == 0)
				return 0;
//				throw new Exception("The arrays could not be empty when computed mean!");
			for (int i = 0; i < size; i++) {
				sum += nums[i];
			}
			mean = Math.round(sum / size);
		} catch (Exception e) {
			logger.error(MathTool.class.getName() + "\t" + e.getMessage());
		}
		return mean;
	}

	public static Integer sd(List<Integer> nums) {
		return MathTool.sd(nums.toArray(new Integer[nums.size()]));
	}

	public static Integer sd(Integer[] nums) {
		int sd = 0;
		try {
			int mean = MathTool.mean(nums);
			int size = nums.length;
			double diff = 0.0f;
			for (int i = 0; i < size; i++) {
				diff += Math.pow(nums[i] - mean, 2);
			}
			sd = (int) Math.round(Math.sqrt(diff / (size - 1)));
		} catch (Exception e) {
			logger.error(MathTool.class.getName() + "\t" + e.getMessage());
		}
		return sd;
	}
	
	public static Integer sd(int[] nums) {
		int sd = 0;
		try {
			int mean = MathTool.mean(nums);
			int size = nums.length;
			double diff = 0.0f;
			for (int i = 0; i < size; i++) {
				diff += Math.pow(nums[i] - mean, 2);
			}
			sd = (int) Math.round(Math.sqrt(diff / (size - 1)));
		} catch (Exception e) {
			logger.error(MathTool.class.getName() + "\t" + e.getMessage());
		}
		return sd;
	} 

	public static Integer avgSd(List<Integer> sds) {
		long sum = 0;
		for (int i = 0; i < sds.size(); i++) {
			sum += Math.pow(sds.get(i), 2);
		}
		int sd = (int) Math.sqrt(sum);
		return sd;
	}
	
	public static int max(List<Integer> values)
	{
		Integer [] arrs = values.toArray(new Integer[values.size()]);
		return max(arrs);
	}
	
	public static int max(Integer[] values)
	{
		Arrays.sort(values);
		return values[values.length - 1];
	}

	public static int max(int typeA, int typeB, int typeC, int typeD) {
		// TODO Auto-generated method stub
		int[] arr = { typeA, typeB, typeC, typeD };
		Arrays.sort(arr);
		return arr[3];
	}

	public static int min(Integer[] values) {
		// TODO Auto-generated method stub
		Arrays.sort(values);
		return values[0];
	}

	public static int min(List<Integer> values) {
		Integer [] arrs = values.toArray(new Integer[values.size()]);
		return min(arrs);
	}

//	public static Map<String, Integer> summary(List<Integer> suptLinks) {
//		int mean = mean(suptLinks);
//		int median = median(suptLinks);
//		int min = min(suptLinks);
//		int max = max(suptLinks);
//		int sd = sd(suptLinks);
//		Map<String, Integer> values = new HashMap<String, Integer>();
//		values.put("MEAN", mean);
//		values.put("MEDIAN", median);
//		values.put("MIN", min);
//		values.put("MAX", max);
//		values.put("SD", sd);
//		return values;
//	}

}
