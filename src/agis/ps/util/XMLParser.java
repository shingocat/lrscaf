/*
*File: agis.ps.util.XMLParser.java
*User: mqin
*Email: mqin@ymail.com
*Date: 2016年1月22日
*/
package agis.ps.util;

import java.io.File;
import java.io.IOException;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

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
			logger.debug(this.getClass().getName() + "\t" + rootElm.getTagName());

			NodeList nodes = rootElm.getChildNodes();
			if (nodes == null || nodes.getLength() == 0) {
				logger.debug(this.getClass().getName() + "\t" + "The XML file do not contain any elements!");
				logger.error(this.getClass().getName() + "\t" + "The XML file do not contain any elements!");
				return null;
			}

			// traverse child input elements
			nodes = rootElm.getElementsByTagName("input");
			if (nodes == null || nodes.getLength() == 0) {
				logger.debug(this.getClass().getName() + "\t" + "The XML file do not contain input elements!");
				logger.error(this.getClass().getName() + "\t" + "The XML file do not contain input elements!");
				return null;
			} else {
				Node input = nodes.item(0);
				nodes = input.getChildNodes();
				for (int i = 0; i < nodes.getLength(); i++) {
					Node node = nodes.item(i);
					if (node.getNodeType() == Node.ELEMENT_NODE) {
						String nodeName = ((Element) node).getNodeName();
						if (nodeName.equalsIgnoreCase("contig")) {
							para.setCntFile(node.getTextContent().trim());
						} else if (nodeName.equalsIgnoreCase("m5")) {
							para.setAlgFile(node.getTextContent().trim());
							para.setType("m");
						} else if (nodeName.equalsIgnoreCase("sam")) {
							para.setAlgFile(node.getTextContent().trim());
							para.setType("s");
						} else if (nodeName.equalsIgnoreCase("bam")) {
							para.setAlgFile(node.getTextContent().trim());
							para.setType("s");
						} else {
							logger.debug(this.getClass().getName() + "\t" + "The XML file contain illeagle elements!");
							logger.error(this.getClass().getName() + "\t" + "The XML file contain illeagle elements!");
							return null;
						}
					}
				}
			}

			// traverse child output elements
			nodes = rootElm.getElementsByTagName("output");
			if (nodes == null || nodes.getLength() == 0) {
				logger.debug(this.getClass().getName() + "\t" + "The XML file do not contain output elements, it will use "
						+ System.getProperty("user.dir") + "!");
				logger.error(this.getClass().getName() + "\t" + "The XML file do not contain output elements, it will use "
						+ System.getProperty("user.dir") + "!");
				para.setOutFolder(System.getProperty("user.dir"));
			} else {
				Node node = nodes.item(0);
				para.setOutFolder(node.getTextContent().trim());
			}

			// traverse child paras elements
			nodes = rootElm.getElementsByTagName("paras");
			if (nodes == null || nodes.getLength() == 0) {
				logger.debug(this.getClass().getName() + "\t" + "The XML file do not contain output elements, it will use default value!");
				logger.error(this.getClass().getName() + "\t" + "The XML file do not contain output elements, it will use default value!");
				para.setMinContLen(2000);
				para.setMinPBLen(5000);
				para.setMinOLLen(1500);
				para.setMinOLRatio(0.8);
				para.setMaxOHLen(500);
				para.setMaxOHRatio(0.1);
				para.setMaxEndLen(500);
				para.setMaxEndRatio(0.1);
				para.setMinSupLinks(3);
				para.setMaxSupLinks(50);
				para.setIdentity(0.8);
				para.setUseOLLink(false);
			} else {
				Node parasNode = nodes.item(0);
				nodes = parasNode.getChildNodes();
				for (int i = 0; i < nodes.getLength(); i++) {
					Node node = nodes.item(i);
					if (node.getNodeType() != Node.ELEMENT_NODE)
						continue;
					String nodeName = node.getNodeName();
					if (nodeName.equalsIgnoreCase("min_contig_length")) {
						para.setMinContLen(Integer.valueOf(node.getTextContent().trim()));
					} else if (nodeName.equalsIgnoreCase("min_pacbio_length")) {
						para.setMinPBLen(Integer.valueOf(node.getTextContent().trim()));
					} else if (nodeName.equalsIgnoreCase("min_overlap_length")) {
						para.setMinOLLen(Integer.valueOf(node.getTextContent().trim()));
					} else if (nodeName.equalsIgnoreCase("min_overlap_ratio")) {
						para.setMinOLRatio(Double.valueOf(node.getTextContent().trim()));
					} else if (nodeName.equalsIgnoreCase("min_supported_links")) {
						para.setMinSupLinks(Integer.valueOf(node.getTextContent().trim()));
					} else if (nodeName.equalsIgnoreCase("max_supported_links")) {
						para.setMaxSupLinks(Integer.valueOf(node.getTextContent().trim()));
					} else if (nodeName.equalsIgnoreCase("max_overhang_length")) {
						para.setMaxOHLen(Integer.valueOf(node.getTextContent().trim()));
					} else if (nodeName.equalsIgnoreCase("max_overhang_ratio")) {
						para.setMaxOHRatio(Double.valueOf(node.getTextContent().trim()));
					} else if (nodeName.equalsIgnoreCase("max_end_length")) {
						para.setMaxEndLen(Integer.valueOf(node.getTextContent().trim()));
					} else if (nodeName.equalsIgnoreCase("max_end_ratio")) {
						para.setMaxEndRatio(Double.valueOf(node.getTextContent().trim()));
					} else if (nodeName.equalsIgnoreCase("identity")) {
						para.setIdentity(Double.valueOf(node.getTextContent().trim()));
					} else if (nodeName.equalsIgnoreCase("use_overlap_link")) {
						String temp = node.getTextContent().trim();
						if(temp.startsWith("t") || temp.startsWith("T"))
							para.setUseOLLink(true);
						else 
							para.setUseOLLink(false);
//						Pattern pat = Pattern.compile("^T", Pattern.CASE_INSENSITIVE);
//						Matcher mat = pat.matcher(node.getTextContent());
//						if(mat.find())
//							para.setUseOLLink(true);
//						else 
//							para.setUseOLLink(false);
					} else {
						logger.debug(this.getClass().getName() + "\t" + "The para element contain illeage item " + nodeName + ". it will omit!");
						logger.info(this.getClass().getName() + "\t" + "The para element contain illeage item " + nodeName + ". it will omit!");
					}
				}
			}

		} catch (ParserConfigurationException e) {
			logger.debug(this.getClass().getName() + "\t" + e.getMessage());
			logger.error(this.getClass().getName() + "\t" + e.getMessage());
		} catch (SAXException e) {
			logger.debug(this.getClass().getName() + "\t" + e.getMessage());
			logger.error(this.getClass().getName() + "\t" + e.getMessage());
		} catch (IOException e) {
			// TODO Auto-generated catch block
			logger.debug(this.getClass().getName() + "\t" + e.getMessage());
			logger.error(this.getClass().getName() + "\t" + e.getMessage());
		}
		return para;
	}
}
