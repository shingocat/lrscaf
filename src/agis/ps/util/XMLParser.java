/*
*File: agis.ps.util.XMLParser.java
*User: mqin
*Email: mqin@ymail.com
*Date: 2016年1月22日
*/
package agis.ps.util;

import java.io.File;
import java.io.IOException;

import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;
import javax.xml.parsers.ParserConfigurationException;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.w3c.dom.Node;
import org.w3c.dom.NodeList;
import org.xml.sax.SAXException;

public class XMLParser {
	private static Logger logger = LoggerFactory.getLogger(XMLParser.class);

	public Parameter parseXML(String file) {
		Parameter para = null;
		try {
			DocumentBuilderFactory factory = DocumentBuilderFactory.newInstance();
			DocumentBuilder builder = factory.newDocumentBuilder();
			Document doc = builder.parse(new File(file));
			para = new Parameter();
			// root element
			Element rootElm = doc.getDocumentElement();
			logger.debug(rootElm.getTagName());

			NodeList nodes = rootElm.getChildNodes();
			if (nodes == null || nodes.getLength() == 0) {
				logger.debug("The XML file do not contain any elements!");
				logger.error("The XML file do not contain any elements!");
				return null;
			}

			// traverse child input elements
			nodes = rootElm.getElementsByTagName("input");
			if (nodes == null || nodes.getLength() == 0) {
				logger.debug("The XML file do not contain input elements!");
				logger.error("The XML file do not contain input elements!");
				return null;
			} else {
				Node input = nodes.item(0);
				nodes = input.getChildNodes();
				for (int i = 0; i < nodes.getLength(); i++) {
					Node node = nodes.item(i);
					if (node.getNodeType() == Node.ELEMENT_NODE) {
						String nodeName = ((Element) node).getNodeName();
						if (nodeName.equalsIgnoreCase("contig")) {
							para.setCntFile(node.getTextContent());
						} else if (nodeName.equalsIgnoreCase("m5")) {
							para.setAlgFile(node.getTextContent());
							para.setType("m");
						} else if (nodeName.equalsIgnoreCase("sam")) {
							para.setAlgFile(node.getTextContent());
							para.setType("s");
						} else if (nodeName.equalsIgnoreCase("bam")) {
							para.setAlgFile(node.getTextContent());
							para.setType("b");
						} else {
							logger.debug("The XML file contain illeagle elements!");
							logger.error("The XML file contain illeagle elements!");
							return null;
						}
					}
				}
			}

			// traverse child output elements
			nodes = rootElm.getElementsByTagName("output");
			if (nodes == null || nodes.getLength() == 0) {
				logger.debug("The XML file do not contain output elements, it will use "
						+ System.getProperty("user.dir") + "!");
				logger.error("The XML file do not contain output elements, it will use "
						+ System.getProperty("user.dir") + "!");
				para.setOutFolder(System.getProperty("user.dir"));
			} else {
				Node node = nodes.item(0);
				para.setOutFolder(node.getTextContent());
			}

			// traverse child paras elements
			nodes = rootElm.getElementsByTagName("paras");
			if (nodes == null || nodes.getLength() == 0) {
				logger.debug("The XML file do not contain output elements, it will use default value!");
				logger.error("The XML file do not contain output elements, it will use default value!");
				para.setMinContLen(1000);
				para.setMaxSupLinks(50);
				para.setMinOLLen(2000);
				para.setMinOLRatio(0.8);
				para.setMinPBLen(5000);
				para.setMinSupLinks(3);
			} else {
				Node parasNode = nodes.item(0);
				nodes = parasNode.getChildNodes();
				for (int i = 0; i < nodes.getLength(); i++) {
					Node node = nodes.item(i);
					if (node.getNodeType() != Node.ELEMENT_NODE)
						continue;
					String nodeName = node.getNodeName();
					if (nodeName.equalsIgnoreCase("min_contig_length")) {
						para.setMinContLen(Integer.valueOf(node.getTextContent()));
					} else if (nodeName.equalsIgnoreCase("min_pacbio_length")) {
						para.setMinPBLen(Integer.valueOf(node.getTextContent()));
					} else if (nodeName.equalsIgnoreCase("min_overlap_length")) {
						para.setMinOLLen(Integer.valueOf(node.getTextContent()));
					} else if (nodeName.equalsIgnoreCase("min_overlap_ratio")) {
						para.setMinOLRatio(Double.valueOf(node.getTextContent()));
					} else if (nodeName.equalsIgnoreCase("min_supported_links")) {
						para.setMinSupLinks(Integer.valueOf(node.getTextContent()));
					} else if (nodeName.equalsIgnoreCase("max_supported_links")) {
						para.setMaxSupLinks(Integer.valueOf(node.getTextContent()));
					} else {
						logger.debug("The para element contain illeage item " + nodeName + ". it will omit!");
						logger.error("The para element contain illeage item " + nodeName + ". it will omit!");
					}
				}
			}

		} catch (ParserConfigurationException e) {
			logger.debug(e.getMessage());
			logger.error(e.getMessage());
		} catch (SAXException e) {
			logger.debug(e.getMessage());
			logger.error(e.getMessage());
		} catch (IOException e) {
			// TODO Auto-generated catch block
			logger.debug(e.getMessage());
			logger.error(e.getMessage());
		}
		return para;
	}
}
