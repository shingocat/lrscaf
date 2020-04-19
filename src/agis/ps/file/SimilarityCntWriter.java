/*
*File: agis.ps.file.SimilarityCntWriter.java
*User: mqin
*Email: mqin@ymail.com
*Date: 2016年12月11日
*/
package agis.ps.file;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.LinkedList;
//import java.util.List;
//import java.util.Map;
//import java.util.Vector;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

//import agis.ps.link.MRecord;
import agis.ps.util.Parameter;

public class SimilarityCntWriter {
	public static Logger logger = LoggerFactory.getLogger(SimilarityCntWriter.class);
	private Parameter paras;
	private File file = null;
	private FileWriter fw = null;
	private BufferedWriter bw = null;

	public SimilarityCntWriter(Parameter paras) {
		this.paras = paras;
		this.init();
	}

	public void init() {
		try {
			file = new File(paras.getOutFolder() + System.getProperty("file.separator") + "simcnts.info");
			if (!file.exists())
				file.createNewFile();
			fw = new FileWriter(file, true);
			bw = new BufferedWriter(fw);
		} catch (IOException e) {
			logger.debug("Error: ", e);
			logger.error(this.getClass().getName() + "\t" + e.getMessage());
		} catch (Exception e) {
			logger.debug("Error: ", e);
			logger.error(this.getClass().getName() + "\t" + e.getMessage());
		}
	}

	public void write(LinkedList<String> sims) {
		try {
			for (String s: sims) {
				if(s.equals(";"))
					bw.newLine();
				else
					bw.write(s + "\t");
//				for(MRecord m : ms)
//				{
//					bw.write(m.gettName() + "\t");
//				}
//				bw.newLine();
			}
			bw.flush();
		} catch (IOException e) {
			logger.debug("Error: ", e);
			logger.error(this.getClass().getName() + "\t" + e.getMessage() + "\t" + e.getClass().getName());
		} finally {
			try {
				if (bw != null)
					bw.close();
			} catch (IOException e) {
				logger.debug("Error: ", e);
				logger.error(this.getClass().getName() + "\t" + e.getMessage() + "\t" + e.getClass().getName());
			}
		}
	}
}
