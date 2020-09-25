package agis.ps.file;
/**
* @author Mao Qin
* @version 2020年9月17日 下午5:49:06
* @Email mqin@outlook.com
* @Description class usage.
* @Copyright All Right Reserved 2020.
*/

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.List;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import agis.ps.path.NodePath;
import agis.ps.util.Parameter;

public class NodePathWriter {

	private static final Logger logger = LoggerFactory.getLogger(NodePathWriter.class);
	private Parameter paras;

	public NodePathWriter(Parameter paras) {
		this.paras = paras;
	}

	public void write(List<NodePath> paths) {
		File out = null;
		BufferedWriter bw = null;
		try {
			out = Paths.get(this.paras.getOutFolder(), "nodePaths.info").toFile();
			// if (out.exists()) {
			// logger.error(filePath + " is exist!");
			// return;
			// }
			if (out.exists()) {
				logger.info("The output file " + out.getCanonicalPath() + " existed. It will overwrite.");
			} else {
				if (!out.createNewFile()) {
					logger.info("The output file" + out.getCanonicalPath() + "could not be created!");
					return;
				}
			}
			bw = Files.newBufferedWriter(out.toPath());
			int count = 0;
			for (NodePath p : paths) {
				bw.write("digraph G" + count + "{\n");
				bw.write(p.toString() + ";\n");
				bw.write("}\n");
				count++;
			}

		} catch (IOException e) {
			logger.debug("Error: ", e);
			logger.error(DotGraphFileWriter.class.getName() + "\t" + e.getMessage());
		} finally {
			if (bw != null)
				try {
					bw.close();
				} catch (IOException e) {
					logger.debug("Error: ", e);
					logger.error(DotGraphFileWriter.class.getName() + "\t" + e.getMessage());
				}
		}
	}
}
