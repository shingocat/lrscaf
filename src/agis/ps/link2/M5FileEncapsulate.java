/*
*File: agis.ps.link2.M5s.java
*User: mqin
*Email: mqin@ymail.com
*Date: 2016年9月23日
*/
package agis.ps.link2;

import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.Map;

import agis.ps.link.M5Record;
import agis.ps.util.Strand;

public class M5FileEncapsulate {
	private int [] qIds;
	private int [] qLengths;
	private int [] qStart;
	private int [] qEnd;
	private int [] tIds;
	private int [] tLengths;
	private int [] tStart;
	private int [] tEnd;
	private byte [] tStrand;
	//private short [] score;
	//private int [] mapQV;
	private int [] matchs;
	private int [] misMatchs;
	private int [] inserts;
	private int [] deletes;
	private int index = 0;
	private int pIdIndex = 0;
	private int cIdIndex = 0;
	private Map<String, Integer> pIdHashForward = new HashMap<String, Integer>();
//	private Map<Integer, String> pIdHashReverse = new HashMap<Integer, String>();
	private Map<String, Integer> cIdHashForward = new LinkedHashMap<String, Integer>();
	private Map<String, Integer> cntLens = new HashMap<String, Integer>();
//	private Map<Integer, String> cIdHashReverse = new HashMap<Integer, String>();
	
	public M5FileEncapsulate()
	{
		
	}
	
	public M5FileEncapsulate(int size)
	{
		this.qIds = new int[size];
		this.qLengths = new int[size];
		this.qStart = new int[size];
		this.qEnd = new int[size];
		this.tIds = new int[size];
		this.tLengths = new int[size];
		this.tStart = new int[size];
		this.tEnd = new int[size];
		this.tStrand = new byte[size];
		this.matchs = new int[size];
		this.misMatchs = new int[size];
		this.inserts = new int[size];
		this.deletes = new int[size];
	}
	
	public M5Record getM5Record(int index)
	{
		if(index < 0 || index >= qIds.length)
			return null;
		M5Record m5 = new M5Record();
		m5.setqName(String.valueOf(qIds[index]));
		m5.setqLength(qLengths[index]);
		m5.setqStart(qStart[index]);
		m5.setqEnd(qEnd[index]);
		m5.settName(String.valueOf(tIds[index]));
		m5.settLength(tLengths[index]);
		m5.settStart(tStart[index]);
		m5.settEnd(tEnd[index]);
		m5.settStrand(tStrand[index] == 0 ? Strand.FORWARD : Strand.REVERSE);
		m5.setNumMatch(matchs[index]);
		m5.setNumMismatch(misMatchs[index]);
		m5.setNumIns(inserts[index]);
		m5.setNumDel(deletes[index]);
		return m5;
	}
	
	// using init array size
	public void addByArraySize(M5Record m5)
	{
		this.addRecord(m5, index);
		index++;
	}
	
	
	// using arrayCopy method
	public void addByArrayCopy(M5Record m5)
	{
		int size = index + 1;
		if(size == 1)
		{
			qIds = new int[size];
			qLengths = new int[size];
			qStart = new int[size];
			qEnd = new int[size];
			tIds = new int[size];
			tLengths = new int[size];
			tStart = new int[size];
			tEnd = new int[size];
			tStrand = new byte[size];
			matchs = new int[size];
			misMatchs = new int[size];
			inserts = new int[size];
			deletes = new int[size];
			this.addRecord(m5, 0);
		} else
		{
			// temp arrays to store data;
			int [] tqIds = new int[size];
			int [] tqLengths = new int[size];
			int [] tqStart = new int[size];
			int [] tqEnd = new int[size];
			int [] ttIds = new int [size];
			int [] ttLengths = new int[size];
			int [] ttStart = new int[size];
			int [] ttEnd = new int[size];
			byte [] ttStrand = new byte[size];
			int [] tmatchs = new int[size];
			int [] tmisMatchs = new int[size];
			int [] tinserts = new int[size];
			int [] tdeletes = new int[size];
			// copy former array into new temp array;
			System.arraycopy(qIds, 0, tqIds, 0, index);
			System.arraycopy(qLengths, 0, tqLengths, 0, index);
			System.arraycopy(qStart, 0, tqStart, 0, index);
			System.arraycopy(qEnd, 0, tqEnd, 0, index);
			System.arraycopy(tIds, 0, ttIds, 0, index);
			System.arraycopy(tLengths, 0, ttLengths, 0, index);
			System.arraycopy(tStart, 0, ttStart, 0, index);
			System.arraycopy(tEnd, 0, ttEnd, 0, index);
			System.arraycopy(tStrand, 0, ttStrand, 0, index);
			System.arraycopy(matchs, 0, tmatchs, 0, index);
			System.arraycopy(misMatchs, 0, tmisMatchs, 0, index);
			System.arraycopy(inserts, 0, tinserts, 0, index);
			System.arraycopy(deletes, 0, tdeletes, 0, index);
			// reassign to origin value;
			qIds = tqIds;
			qLengths = tqLengths;
			qStart = tqStart;
			qEnd = tqEnd;
			tIds = ttIds;
			tLengths = ttLengths;
			tStart = ttStart;
			tEnd = ttEnd;
			tStrand = ttStrand;
			matchs = tmatchs;
			misMatchs = tmisMatchs;
			inserts = tinserts;
			deletes = tdeletes;
			this.addRecord(m5, size - 1);
		}
		index++;
	}
	
	public void addRecord(M5Record m5, int index)
	{
		// for the last element to insert into array
		String pId = m5.getqName();
		int pIdInt = 0;
		String cId = m5.gettName();
		int cIdInt = 0;
		if(pIdHashForward.containsKey(pId))
		{
			pIdInt = pIdHashForward.get(pId);
		} else
		{
			pIdInt = pIdIndex;
			pIdHashForward.put(pId, pIdIndex);
			pIdIndex++;
		}
		if(cIdHashForward.containsKey(cId))
		{
			cIdInt = cIdHashForward.get(cId);
		} else
		{
			cIdInt = cIdIndex;
			cIdHashForward.put(cId, cIdIndex);
			cIdIndex++;
		}
		qIds[index] = pIdInt;
		qLengths[index] = m5.getqLength();
		qStart[index] = m5.getqStart();
		qEnd[index] = m5.getqEnd();
		tIds[index] = cIdInt;
		tLengths[index] = m5.gettLength();
		tStart[index] = m5.gettStart();
		tEnd[index] = m5.gettEnd();
		if(m5.gettStrand().equals(Strand.FORWARD))
			tStrand[index] = 0;
		else
			tStrand[index] = 1;
		matchs[index] = m5.getNumMatch();
		misMatchs[index] = m5.getNumMismatch();
		inserts[index] = m5.getNumIns();
		deletes[index] = m5.getNumDel();
		if(!cntLens.containsKey(cId))
			cntLens.put(cId, m5.gettLength());
		// next array element index;
//		index++;
	}
	
	public String getCntId(int index)
	{
		if(cIdHashForward.containsValue(index))
		{
			return null;
		} else
		{
			return null;
		}
	}
	
	public CntFileEncapsulate getCntFileEncapsulate()
	{
		CntFileEncapsulate cntFile = new CntFileEncapsulate();
		Iterator<String> it = cIdHashForward.keySet().iterator();
		while(it.hasNext())
		{
			String oid = it.next();
			Integer nid = this.cIdHashForward.get(oid);
			int len = cntLens.get(oid);
			cntFile.addLength(nid, oid, len);
		}
		return cntFile;
	}
}


