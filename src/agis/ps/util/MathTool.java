/*
*File: agis.ps.util.MathTool.java
*User: mqin
*Email: mqin@ymail.com
*Date: 2016年1月13日
*/
package agis.ps.util;

import java.util.Arrays;
import java.util.List;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

public class MathTool {
	private static Logger logger = LoggerFactory.getLogger(MathTool.class);

	public static Integer sum(List<Integer> nums) {
		int sum = 0;
		for (int i = 0; i < nums.size(); i++) {
			sum += nums.get(i);
		}
		return sum;
	}

	public static Integer mean(List<Integer> nums) {
		return MathTool.mean(nums.toArray(new Integer[nums.size()]));
	}

	public static Integer mean(Integer[] nums) {
		int mean = 0;
		try {
			int sum = 0;
			int size = nums.length;
			if (size == 0)
				throw new Exception("The arrays could not be empty when computed mean!");
			for (int i = 0; i < size; i++) {
				sum += nums[i];
			}
			mean = Math.round(sum / size);
		} catch (Exception e) {
			logger.debug(MathTool.class.getName() + "\t" + e.getMessage());
			logger.info(MathTool.class.getName() + "\t" + e.getMessage());
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
			logger.debug(MathTool.class.getName() + "\t" + e.getMessage());
			logger.info(MathTool.class.getName() + "\t" + e.getMessage());
		}
		return sd;
	}

	public static Integer avgSd(List<Integer> sds) {
		int sum = 0;
		for (int i = 0; i < sds.size(); i++) {
			sum += Math.pow(sds.get(i), 2);
		}
		int sd = (int) Math.sqrt(sum);
		return sd;
	}

	public static int max(int typeA, int typeB, int typeC, int typeD) {
		// TODO Auto-generated method stub
		int[] arr = { typeA, typeB, typeC, typeD };
		Arrays.sort(arr);
		return arr[3];
	}

	public static int min(Integer[] values) {
		// TODO Auto-generated method stub

		return 0;
	}

	public static int min(List<Integer> values) {
		return 0;
	}

}
