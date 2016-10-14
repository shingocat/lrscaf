/*
*File: agis.ps.file.OutputFolderBuilder.java
*User: mqin
*Email: mqin@ymail.com
*Date: 2016年7月19日
*/
package agis.ps.file;

import java.io.File;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import agis.ps.util.Parameter;

public class OutputFolderBuilder {

	private static Logger logger = LoggerFactory.getLogger(OutputFolderBuilder.class);
	private Parameter paras;

	public OutputFolderBuilder(Parameter paras) {
		this.paras = paras;
	}

	public boolean building() {
		long start = System.currentTimeMillis();
		String path = paras.getOutFolder();
		boolean isValid = false;
		if (path == null || path.length() == 0) {
			logger.error(this.getClass().getName() + "\t" + "The output path was not setted!");
			return isValid;
		}
		try {
			File output = new File(path);
			if (output.exists()) {
				logger.info(this.getClass().getName() + "\t" + "The output folder was exist!");
				logger.info(this.getClass().getName() + "\t" + "It will delete all file under this folder!");
				isValid = this.deleteDir(output);
			} else {
				isValid = output.mkdirs();
				if (isValid)
					logger.info(this.getClass().getName() + "\t" + "Build output folder successfully!");
				else
					logger.info(this.getClass().getName() + "\t" + "Build output folder failed!");
			}
		} catch (SecurityException e) {
			logger.error(this.getClass().getName() + e.getMessage() + "\t" + e.getClass().getName());
		} catch (Exception e) {
			logger.error(this.getClass().getName() + e.getMessage() + "\t" + e.getClass().getName());
		}
		long end = System.currentTimeMillis();
		logger.info("Building output folder, erase time: " + (end - start) + " ms");
		return isValid;
	}

	private boolean deleteDir(File dir) {
		boolean success = false;
		if (dir.isDirectory()) {
			String[] children = dir.list();
			if(children.length == 0)
				success = true;
			for(String s : children)
			{
				File file = new File(dir,s);
				if(file.isFile())
				{
					success = file.delete();
					if(!success)
						break;
				}
			}
//			for (int i = 0; i < children.length; i++) {
//				boolean success = deleteDir(new File(dir, children[i]));
//				if (!success) {
//					return false;
//				}
//			}
		}
//		return dir.delete();
		return success;
	}
}
