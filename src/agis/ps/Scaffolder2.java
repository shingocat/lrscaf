/*
*File: agis.ps.Scaffolder2.java
*User: mqin
*Email: mqin@ymail.com
*Date: 2016年7月17日
*/
package agis.ps;

import java.io.File;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import agis.ps.util.Parameter;

public class Scaffolder2 {
	private static Logger logger = LoggerFactory.getLogger(Scaffolder.class);
	private Parameter paras;
	
	public Scaffolder2(Parameter paras)
	{
		this.paras = paras;
	}
	
	public void scaffolding()
	{
		long start = System.currentTimeMillis();
		logger.info("starting...");
		try
		{
			// building output folder
			if(!buildOutputFolder())
				return;
			
		} catch(Exception e)
		{
			logger.error(this.getClass().getName() + "\t" + e.getMessage());
		}
		logger.info("Ending...");
		long end = System.currentTimeMillis();
		logger.info("Scaffolding erase time: " + (end - start)/1000 + " s");
	}
	
	// building output folder 
	private boolean buildOutputFolder()
	{
		long start = System.currentTimeMillis();
		String path = paras.getOutFolder();
		boolean isValid = false;
		if (path == null || path.length() == 0) {
			logger.error(this.getClass().getName() + "The output path was not setted!");
			logger.debug(this.getClass().getName() + "The output path was not setted!");
			return isValid;
		}
		try {
			File output = new File(path);
			if (output.exists()) {
				logger.info(this.getClass().getName() + "The output folder was exist!");
				logger.debug(this.getClass().getName() + "The output folder was exist!");
				isValid = true;
			} else {
				isValid = output.mkdirs();
				if(isValid)
				{
					logger.info(this.getClass().getName() + "\t" + "Build output folder successfully!");
					logger.debug(this.getClass().getName() + "\t" + "Build output folder successfully!");
				} else
				{
					logger.info(this.getClass().getName() + "\t" + "Build output folder failed!");
					logger.debug(this.getClass().getName() + "\t" + "Build output folder failed!");
				}
			}
		} catch (SecurityException e) {
			logger.error(this.getClass().getName() + e.getMessage() + "\t" + e.getClass().getName());
			logger.debug(this.getClass().getName() + e.getMessage() + "\t" + e.getClass().getName());
		} catch (Exception e) {
			logger.error(this.getClass().getName() + e.getMessage() + "\t" + e.getClass().getName());
			logger.debug(this.getClass().getName() + e.getMessage() + "\t" + e.getClass().getName());
		}
		long end = System.currentTimeMillis();
		logger.info("Building output folder, erase time: " + (end - start)/1000 + " s");
		return isValid;
	}
	
}


