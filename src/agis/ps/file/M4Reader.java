/*
*File: agis.ps.file.M4Reader.java
*User: mqin
*Email: mqin@ymail.com
*Date: 2016年5月30日
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

import agis.ps.link.M4Record;
import agis.ps.link.MRecord;
import agis.ps.util.Parameter;
import agis.ps.util.Strand;

public class M4Reader {

	private final static Logger logger = LoggerFactory.getLogger(M4Reader.class);
	private String path;
	private List<MRecord> m4List;
	private int minPBLen = 5000; // minimum pacbio read length
	private int minCNTLen = 3000; // minimum contig length;
	private double identity = 0.8;
	
	public M4Reader(Parameter paras)
	{
		this.path = paras.getAlgFile();
		this.minPBLen = paras.getMinPBLen();
		this.minCNTLen = paras.getMinContLen();
		this.identity = paras.getIdentity();
	}
	
	
	public List<MRecord> read()
	{
		long start = System.currentTimeMillis();
		FileReader fr = null;
		BufferedReader br = null;
		int count = 0;
		try
		{
			File mFile = new File(path);
			if(!mFile.exists())
			{
				logger.debug(this.getClass().getName() + "The m4 file do not exist!");
				logger.error(this.getClass().getName() + "The m4 file do not exist!");
				return null;
			}
			
			fr = new FileReader(mFile);
			br = new BufferedReader(fr);
			String line = null;
			String [] arrs = null;
			m4List = new Vector<MRecord>();
			while((line = br.readLine()) != null)
			{
				line = line.trim();
				line = line.replaceAll(System.getProperty("line.separator"), "");
				arrs = line.split("\\s+");
				if(arrs[0].equalsIgnoreCase("qName") && arrs[1].equalsIgnoreCase("tName"))
					continue;
				count++;
				// if less than minimum pacbio length
				if(Integer.valueOf(arrs[7]) < minPBLen)
					continue;
				// if less tan minimum contig length
				if(Integer.valueOf(arrs[11]) < minCNTLen)
					continue;
				// if the identity less than specified value;
				double value = Double.valueOf(arrs[3]);
				if(value < identity)
					continue;
				M4Record m4 = new M4Record();
				m4.setqName(arrs[0]);
				m4.setqLength(Integer.valueOf(arrs[7]));
				m4.setqStart(Integer.valueOf(arrs[5]));
				m4.setqEnd(Integer.valueOf(arrs[6]));
				m4.setqStrand(arrs[4].equals("0") ? Strand.FORWARD : Strand.REVERSE);
				m4.settName(arrs[1]);
				m4.settLength(Integer.valueOf(arrs[11]));
				// Reverse strand alignments for m 4 are reported on the reverse strand (3'->5')
				if(arrs[8].equals("1"))
				{
					m4.settStart(m4.gettLength() - Integer.valueOf(arrs[10]));
					m4.settEnd(m4.gettLength() - Integer.valueOf(arrs[9]));
				} else
				{
					m4.settStart(Integer.valueOf(arrs[9]));
					m4.settEnd(Integer.valueOf(arrs[10]));
				}
				m4.settStrand(arrs[8].equals("0") ? Strand.FORWARD : Strand.REVERSE);
				m4.setScore(Integer.valueOf(arrs[2]));
				m4.setMapQV(Integer.valueOf(arrs[12]));
				// setting identity
				m4.setIdentity(value);
				m4List.add(m4);
			}
			br.close();
		} catch(ArrayIndexOutOfBoundsException e)
		{
			logger.debug(this.getClass().getName() + e.getMessage() + "\t" + e.getClass().getName());
			logger.error(this.getClass().getName() + e.getMessage() + "\t" + e.getClass().getName());
			m4List = null;
		} catch(FileNotFoundException e)
		{
			logger.debug(this.getClass().getName() + e.getMessage() + "\t" + e.getClass().getName());
			logger.error(this.getClass().getName() + e.getMessage() + "\t" + e.getClass().getName());
			m4List = null;
		} catch(IOException e)
		{
			logger.debug(this.getClass().getName() + e.getMessage() + "\t" + e.getClass().getName());
			logger.error(this.getClass().getName() + e.getMessage() + "\t" + e.getClass().getName());
			m4List = null;
		} catch(Exception e){
			logger.debug(this.getClass().getName() + e.getMessage() + "\t" + e.getClass().getName());
			logger.error(this.getClass().getName() + e.getMessage() + "\t" + e.getClass().getName());
			m4List = null;
		} finally
		{
			try
			{
				if( br != null)
					br.close();
			} catch(IOException e)
			{
				logger.debug(this.getClass().getName() + e.getMessage() + "\t" + e.getClass().getName());
				logger.error(this.getClass().getName() + e.getMessage() + "\t" + e.getClass().getName());
				m4List = null;
			}
		}
		long end = System.currentTimeMillis();
		logger.info(this.getClass().getName() + "\tAligned records: " + count);
		logger.info(this.getClass().getName() + "\tRecords after filtered: " + m4List.size());
		logger.info("Reading M4 Record, erase time: " + (end - start) + " ms");
		return m4List;
	}
}


