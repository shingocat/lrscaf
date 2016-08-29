/*
*File: agis.ps.file.TriadLinkWriter.java
*User: mqin
*Email: mqin@ymail.com
*Date: 2016年5月13日
*/
package agis.ps.file;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Arrays;
import java.util.List;
import java.util.Map;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import agis.ps.link.Edge;
import agis.ps.link.TriadLink;
import agis.ps.seqs.Contig;
import agis.ps.util.Parameter;

public class TriadLinkWriter {
	public static Logger logger = LoggerFactory.getLogger(TriadLinkWriter.class);
	private Parameter paras;
	private List<TriadLink> triads;
	private File file = null;
	private FileWriter fw = null;
	private BufferedWriter bw = null;

	public TriadLinkWriter(Parameter paras) {
		this.paras = paras;
	}

	public TriadLinkWriter(Parameter paras, List<TriadLink> triads) {
		this.paras = paras;
		this.triads = triads;
	}

	public void init() {
		try {
			file = new File(paras.getOutFolder() + System.getProperty("file.separator") + "triadlinks.info");
			if (!file.exists())
				file.createNewFile();
			fw = new FileWriter(file, true);
			bw = new BufferedWriter(fw);
		} catch (IOException e) {
			logger.error(this.getClass().getName() + "\t" + e.getMessage());
		} catch (Exception e) {
			logger.error(this.getClass().getName() + "\t" + e.getMessage());
		}
	}

	public void write() {
		String outFolder = paras.getOutFolder();
		String fileName = outFolder + System.getProperty("file.separator") + "triad_link.info";
		File file = null;
		FileWriter fw = null;
		BufferedWriter bw = null;
		try {
			file = new File(fileName);
			if (file.exists()) {
				logger.debug(this.getClass().getName() + " The output file " + fileName
						+ " is exist! It will not be overwrited!");
				logger.info(this.getClass().getName() + " The output file " + fileName
						+ " is exist! It will not be overwrited!");
				return;
			}
			if (!file.createNewFile()) {
				logger.debug(this.getClass().getName() + " The output file " + fileName + " could not be created!");
				logger.info(this.getClass().getName() + " The output file " + fileName + " could not be created!");
				return;
			}
			fw = new FileWriter(file);
			bw = new BufferedWriter(fw);
			for (TriadLink tl : triads) {
				Contig pre = tl.getPrevious();
				Contig mid = tl.getMiddle();
				Contig lst = tl.getLast();
				bw.write(pre.getID() + "," + mid.getID() + "," + lst.getID() + "," + tl.getSupLinks());
				// bw.write(pre.getID() + "(length=" + pre.getLength() + ")," +
				// mid.getID() + "(length=" + mid.getLength() + ")," +
				// lst.getID() + "(length=" + lst.getLength() + ")," +
				// "supLinks=" + tl.getSupLinks());
				bw.newLine();
			}
			bw.flush();
		} catch (IOException e) {
			logger.error(this.getClass().getName() + "\t" + e.getMessage() + "\t" + e.getClass().getName());
		} finally {
			try {
				if (bw != null)
					bw.close();
			} catch (IOException e) {
				logger.error(this.getClass().getName() + "\t" + e.getMessage() + "\t" + e.getClass().getName());
			}
		}
	}

	public void write2(List<TriadLink> triads) {
		try {
			for (TriadLink tl : triads) {
				Contig pre = tl.getPrevious();
				Contig mid = tl.getMiddle();
				Contig lst = tl.getLast();
				String line = null;
				if (mid == null)
					line = pre.getID() + ",-" + "," + lst.getID() + "," + tl.getSupLinks();
				else
					line = pre.getID() + "," + mid.getID() + "," + lst.getID() + "," + tl.getSupLinks();
				bw.write(line);
				bw.newLine();
			}
			bw.flush();
		} catch (IOException e) {
			logger.error(this.getClass().getName() + "\t" + e.getMessage() + "\t" + e.getClass().getName());
		}
	}

	public void write2(TriadLink triad) {
		try {

			Contig pre = triad.getPrevious();
			Contig mid = triad.getMiddle();
			Contig lst = triad.getLast();
			String line = null;
			if (mid == null)
				line = pre.getID() + ",-" + "," + lst.getID() + "," + triad.getSupLinks();
			else
				line = pre.getID() + "," + mid.getID() + "," + lst.getID() + "," + triad.getSupLinks();
			bw.write(line);
			bw.newLine();
			bw.flush();
		} catch (IOException e) {
			logger.error(this.getClass().getName() + "\t" + e.getMessage() + "\t" + e.getClass().getName());
		}
	}
	
	public void write4Edges(List<Edge> edges)
	{
		try
		{
			if(edges != null && !(edges.isEmpty()))
			{
				Edge e = edges.get(0);
				Contig origin = e.getOrigin();
				Contig terminus = e.getTerminus();
				int support = e.getLinkNum();
				String line = origin.getID() + ",-," + terminus.getID() + "," + support;
				bw.write(line);
				bw.newLine();
				bw.flush();
			}

		} catch(IOException e)
		{
			logger.error(this.getClass().getName() + "\t" + e.getMessage() + "\t" + e.getClass().getName());
		} catch(Exception e)
		{
			logger.error(this.getClass().getName() + "\t" + e.getMessage() + "\t" + e.getClass().getName());
		}
	}

	public void close() {
		try {
			if (bw != null)
				bw.close();
		} catch (IOException e) {
			logger.error(this.getClass().getName() + "\t" + e.getMessage());
		}
	}
}
