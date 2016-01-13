/*
*File: agis.ps.util.DotGraphWriter.java
*User: mqin
*Email: mqin@ymail.com
*Date: 2015年12月28日
*/
package agis.ps.util;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import agis.ps.Edge;
import agis.ps.Path;
import agis.ps.SimplePath;
import agis.ps.link.PBLink;

public class DotGraphFileWriter {

	public static final Logger logger = LoggerFactory.getLogger(DotGraphFileWriter.class);

	private String filePath;
	private List<SimplePath> paths;
	private List<PBLink> links;

	public DotGraphFileWriter(String filePath, List<SimplePath> paths) {
		this.filePath = filePath;
		this.paths = paths;
	}
	
//	public DotGraphFileWriter(String filePath, List<PBLink> links)
//	{
//		this.filePath = filePath;
//		this.links = links;
//	}
//	
	public static void writePBLink(String filePath, List<PBLink> links)
	{
		if(filePath == null)
			return;
		File out = null;
		FileWriter fw = null;
		BufferedWriter bw = null;
		try {
			out = new File(filePath);
			if (out.exists()) {
				logger.debug("The output file of dot graph is exist!");
				logger.info("The output file of dot graph is exist!");
				return;
			}
			out.createNewFile();
			fw = new FileWriter(out);
			bw = new BufferedWriter(fw);
			bw.write("digraph G{\n");
			for(PBLink pb : links)
			{
				String color = "red";
				if(pb.getoStrand().equals(Strand.FORWARD) && pb.gettStrand().equals(Strand.REVERSE))
					color = "green";
				else if(pb.getoStrand().equals(Strand.REVERSE) && pb.gettStrand().equals(Strand.REVERSE))
					color = "blue";
				else if(pb.getoStrand().equals(Strand.REVERSE) && pb.gettStrand().equals(Strand.FORWARD))
					color = "yellow";
				bw.write(pb.getOrigin().getID() + " -> " + pb.getTerminus().getID() + " [label=\"" +
						 pb.getoStrand().toString() + pb.getoStartLoc() + ":" + pb.gettStrand().toString() +
						 pb.gettStartLoc() + "\",color=" + color + "];\n");
			}
			bw.write("}");
		} catch (IOException e) {
			logger.debug(e.getMessage());
			logger.error(e.getMessage());
		} finally
		{
			if(bw != null)
				try {
					bw.close();
				} catch (IOException e) {
					logger.debug(e.getMessage());
					logger.error(e.getMessage());;
				}
		}

	}
	
	public static void writeEdge(String filePath, List<Edge> edges)
	{
		if(filePath == null)
			return;
		File out = null;
		FileWriter fw = null;
		BufferedWriter bw = null;
		try {
			out = new File(filePath);
			if (out.exists()) {
				logger.debug("The output file of dot graph is exist!");
				logger.info("The output file of dot graph is exist!");
				return;
			}
			out.createNewFile();
			fw = new FileWriter(out);
			bw = new BufferedWriter(fw);
			bw.write("digraph G{\n");
			for(Edge e : edges)
			{
				bw.write(e.getOrigin().getID() + " -> " + e.getTerminus().getID() + " [label=\"" +
						e.getLinkNum()
						+ "\",color= blue];\n");
			}
			bw.write("}");
		} catch (IOException e) {
			logger.debug(e.getMessage());
			logger.error(e.getMessage());
		} finally
		{
			if(bw != null)
				try {
					bw.close();
				} catch (IOException e) {
					logger.debug(e.getMessage());
					logger.error(e.getMessage());;
				}
		}

	}
	
	public void write() {
		if(paths == null)
			return;
		
		File out = null;
		FileWriter fw = null;
		BufferedWriter bw = null;
		try {
			out = new File(filePath);
			if (out.exists()) {
				logger.debug("The output file of dot graph is exist!");
				logger.info("The output file of dot graph is exist!");
				return;
			}
			out.createNewFile();
			fw = new FileWriter(out);
			bw = new BufferedWriter(fw);
			bw.write("digraph G{\n");
			//for the original paths code
			for(SimplePath sp : paths)
			{
				bw.write(sp.getStart() + " -> " + sp.getEnd() + " [label=\"" +
						sp.getLabel() + "\",color=" + sp.getColor().toString() + "];\n");
			}
			bw.write("}");
		} catch (IOException e) {
			logger.debug(e.getMessage());
			logger.error(e.getMessage());
		} finally
		{
			if(bw != null)
				try {
					bw.close();
				} catch (IOException e) {
					logger.debug(e.getMessage());
					logger.error(e.getMessage());;
				}
		}

	}

}
