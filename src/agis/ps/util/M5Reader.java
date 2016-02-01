/*
*File: agis.ps.util.M5Reader.java
*User: mqin
*Email: mqin@ymail.com
*Date: 2015年12月28日
*/
package agis.ps.util;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Vector;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import agis.ps.M5Record;

public class M5Reader {
	
	private final static Logger logger = LoggerFactory.getLogger(M5Reader.class);
	
	private String path;
	private List<M5Record> m5List;
	
	public M5Reader(String path)
	{
		this.path = path;
	}
	
	public List<M5Record> read()
	{
		FileReader fr = null;
		BufferedReader br = null;
		try
		{
			File m5File = new File(path);
			if(!m5File.exists())
			{
				logger.debug("The m5 file do not exist!");
				logger.error("The m5 file do not exist!");
				return null;
			}
			
			fr = new FileReader(m5File);
			br = new BufferedReader(fr);
			String line = null;
			String [] arrs = null;
			m5List = new Vector<M5Record>();
			while((line = br.readLine()) != null)
			{
				arrs = line.split("\\s+");
				// if the first line is header
				if(arrs[0].equalsIgnoreCase("qName") && arrs[1].equalsIgnoreCase("qLength"))
					continue;
				M5Record m5 = new M5Record();
				m5.setqName(arrs[0]);
				m5.setqLength(Integer.valueOf(arrs[1]));
				m5.setqStart(Integer.valueOf(arrs[2]));
				m5.setqEnd(Integer.valueOf(arrs[3]));
				m5.setqStrand(arrs[4].equals("+") ? Strand.FORWARD : Strand.REVERSE);
				m5.settName(arrs[5]);
				m5.settLength(Integer.valueOf(arrs[6]));
				m5.settStart(Integer.valueOf(arrs[7]));
				m5.settEnd(Integer.valueOf(arrs[8]));
				m5.settStrand(arrs[9].equals("+") ? Strand.FORWARD : Strand.REVERSE);
				m5.setScore(Integer.valueOf(arrs[10]));
				m5.setNumMatch(Integer.valueOf(arrs[11]));
				m5.setNumMismatch(Integer.valueOf(arrs[12]));
				m5.setNumIns(Integer.valueOf(arrs[13]));
				m5.setNumDel(Integer.valueOf(arrs[14]));
				m5.setMapQV(Integer.valueOf(arrs[15]));
				m5.setqAlignedSeq(arrs[16]);
				m5.setMatchPattern(arrs[17]);
				m5.settAlignedSeq(arrs[18]);
				m5List.add(m5);
			}
			br.close();
			fr.close();
		} catch(ArrayIndexOutOfBoundsException e)
		{
			logger.debug(e.getMessage());
			logger.error(e.getMessage());
		} catch(FileNotFoundException e)
		{
			logger.debug(e.getMessage());
			logger.error(e.getMessage());
		} catch(IOException e)
		{
			logger.debug(e.getMessage());
			logger.error(e.getMessage());
		} finally
		{
			if(br != null)
			{
				try {
					br.close();
				} catch (IOException e) {
					logger.debug(e.getMessage());
					logger.error(e.getMessage());
				}
			} 
			if(fr != null)
			{
				try {
					fr.close();
				} catch (IOException e) {
					logger.debug(e.getMessage());
					logger.error(e.getMessage());
				}
			}
		}
		return m5List;
	}
}


