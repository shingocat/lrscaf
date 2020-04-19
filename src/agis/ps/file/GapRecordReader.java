/*
*File: agis.ps.file.GapRecordReader.java
*User: mqin
*Email: mqin@ymail.com
*Date: 2016年6月2日
*/
package agis.ps.file;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.List;
import java.util.Vector;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

//import agis.ps.link.M4Record;
//import agis.ps.link.MRecord;
import agis.ps.seqs.PBGapSeq;
import agis.ps.util.GapRecord;
import agis.ps.util.Parameter;
import agis.ps.util.Strand;

public class GapRecordReader {
	private final static Logger logger = LoggerFactory.getLogger(GapRecordWriter.class);
//	private Parameter paras;
	private String path;
	private List<GapRecord> gaps;
	
	public GapRecordReader(String path)
	{
		this.path = path;
	}
	
	public GapRecordReader(Parameter paras)
	{
//		this.paras = paras;
		this.path = paras.getOutFolder() + System.getProperty("file.separator") + "gap_record.info";
	}
	
	
	public List<GapRecord> read()
	{
		FileReader fr = null;
		BufferedReader br = null;
//		int count = 0;
		try
		{
			File mFile = new File(path);
			if(!mFile.exists())
			{
				logger.debug(this.getClass().getName() + "The GapRecord file do not exist!");
				logger.error(this.getClass().getName() + "The GapRecord file do not exist!");
				return null;
			}
			
			fr = new FileReader(mFile);
			br = new BufferedReader(fr);
			String line = null;
			String [] arrs = null;
			gaps = new Vector<GapRecord>(50);
			GapRecord gr = null;
			while((line = br.readLine()) != null)
			{
				line = line.trim();
				line = line.replaceAll(System.getProperty("line.separator"), "");
				arrs = line.split("\\s+");
				if(arrs[0].equalsIgnoreCase(">"))
				{
					if(gr != null)
					{
						gaps.add(gr);
						gr = null;
					}
					gr = new GapRecord();
					gr.setStart(arrs[1]);
					gr.setEnd(arrs[2]);
				} else
				{
					PBGapSeq seq = new PBGapSeq();
					seq.setId(arrs[0]);
					seq.setStart(Integer.valueOf(arrs[1]));
					seq.setEnd(Integer.valueOf(arrs[2]));
					seq.setStrand(arrs[3].equalsIgnoreCase("+")?Strand.FORWARD:Strand.REVERSE);
					gr.addSeq(seq);
				}
			}
			if(gr != null)
			{
				gaps.add(gr);
				gr = null;
			}
			br.close();
		} catch(ArrayIndexOutOfBoundsException e)
		{
			logger.debug("Error: ", e);
			logger.error(this.getClass().getName() + e.getMessage() + "\t" + e.getClass().getName());
		} catch(FileNotFoundException e)
		{
			logger.debug("Error: ", e);
			logger.error(this.getClass().getName() + e.getMessage() + "\t" + e.getClass().getName());
		} catch(IOException e)
		{
			logger.debug("Error: ", e);
			logger.error(this.getClass().getName() + e.getMessage() + "\t" + e.getClass().getName());
		} catch(Exception e){
			logger.debug("Error: ", e);
			logger.error(this.getClass().getName() + e.getMessage() + "\t" + e.getClass().getName());
		} finally
		{
			try
			{
				if( br != null)
					br.close();
			} catch(IOException e)
			{
				logger.debug("Error: ", e);
				logger.error(this.getClass().getName() + e.getMessage() + "\t" + e.getClass().getName());
			}
		}
		logger.debug(this.getClass().getName() + "\tGapRecord: " + gaps.size());
		return gaps;
	}
}


