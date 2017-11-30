/*
*File: agis.ps.file.N50Writer.java
*User: mqin
*Email: mqin@ymail.com
*Date: 2017年2月8日
*/
package agis.ps.file;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.List;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import agis.ps.util.MathTool;
import agis.ps.util.Parameter;

public class N50Writer {
	private static Logger logger = LoggerFactory.getLogger(N50Writer.class);
	private Parameter paras;
	private List<Integer>lens;
	private String prefix;
	
	
	public N50Writer(Parameter paras, String prefix, List<Integer> lens)
	{
		this.paras = paras;
		this.prefix = prefix;
		this.lens = lens;
	}
	
	public void write()
	{
		if(lens == null)
			throw new IllegalArgumentException(this.getClass() + "\t:The list of sequence length recould could not be null!");
		String outFolder = paras.getOutFolder();
		String fileName = outFolder + System.getProperty("file.separator") + prefix + "_summary.info";
		File file = null; 
		FileWriter fw = null;
		BufferedWriter bw = null;
		try
		{
			file = new File(fileName);
			if (file.exists()) {
				logger.info("The output file" + fileName + "is exist! It will not be overwrited!");
				return;
			}
			if(!file.createNewFile())
			{
				logger.info(this.getClass().getName() + "The output file" + fileName + "could not be created!");
				return;
			}
			fw = new FileWriter(file);
			bw = new BufferedWriter(fw);
			long total = MathTool.sum(lens);
			int mean = MathTool.mean(lens);
			int num = lens.size();
			int median = (int) MathTool.median(lens, false);
			bw.write("Total Length:\t" + total);
			bw.newLine();
			bw.write("Sequences No.:\t" + num);
			bw.newLine();
			bw.write("Mean Length:\t" + mean);
			bw.newLine();
			bw.write("Median:\t" + median);
			bw.newLine();
//			Arrays.sort(lens);
			Collections.sort(lens);
			double index = 0.0;
			long values = 0L;
			long threshold = (long) (index * total);
			int count = 0;
			for(int i = num - 1; i >= 0; i--)
			{
				values += lens.get(i);
				while(values >= threshold)
				{
					bw.write("N" + count + "0:\t" + lens.get(i));
					bw.newLine();
					index += 0.1;
					count++;
					threshold = (long) (index * total);
				}
			}
			bw.flush();
		} catch(IOException e)
		{
			logger.error(this.getClass().getName() + "\t" + e.getMessage() + "\t" + e.getClass().getName());
		} finally
		{
			try{
				if(bw != null)
					bw.close();
			} catch(IOException e)
			{
				logger.error(this.getClass().getName() + "\t" + e.getMessage() + "\t" + e.getClass().getName());
			}
		}
	}
}


