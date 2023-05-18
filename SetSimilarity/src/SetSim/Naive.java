package SetSim;

import java.math.BigInteger;
import java.text.BreakIterator;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;


import SetSim.KDTree.EncLeaf;
import SetSim.Trapdoor.FilterTrapdoor;
import SetSim.Trapdoor.RefinementTrapdoor;


public class Naive {
	
	
	// security parameters
	
	public static SHE paramters = MainFinal.paramters;
	public static BigInteger cipherZero1 = MainFinal.cipherZero1;
	public static BigInteger cipherZero2 = MainFinal.cipherZero2; 
	public static BigInteger cipherMinusOne = MainFinal.cipherMinusOne;
	
	public static HashMap<Integer, ArrayList<BigInteger>> encDataSetNaive(int n, int l) throws Exception {
		
		HashMap<Integer, double[]> dataList = Data.ReadDataJeterWithoutDup(n, l);
		
		HashMap<Integer, ArrayList<Integer>> dataHashMap = DataOutsource.vectorMapToSetList(dataList);
		
		HashMap<Integer, ArrayList<BigInteger>> encDataset = new HashMap<Integer, ArrayList<BigInteger>>();
		
		
		// encrypt dataset
		for (Iterator iterator = dataHashMap.keySet().iterator(); iterator.hasNext();) {
			
			int key = (int) iterator.next();
			ArrayList<Integer> set = dataHashMap.get(key);
			
			ArrayList<BigInteger> encRecord = SHE.encIntegerList(set, paramters);
			encDataset.put(key, encRecord);
			
		}
		
		return encDataset;
		
	}
	
	
	public static HashMap<Integer, ArrayList<BigInteger>> queryNaive(HashMap<Integer, ArrayList<BigInteger>> encDataset, RefinementTrapdoor refinementTrapdoor, 
			HashMap<Integer, ArrayList<BigInteger>> rList, HashMap<Integer, ArrayList<BigInteger>> rAddmList) {
		
		
		for (Iterator iterator = encDataset.keySet().iterator(); iterator.hasNext();) {
			
			int key = (int) iterator.next();
			ArrayList<BigInteger> encX = encDataset.get(key);
			
			if(KDTree.leafEval(encX, refinementTrapdoor, paramters)) {
				
				ArrayList<BigInteger> r = KDTree.genRandomVector(encX.size(), paramters);
				ArrayList<BigInteger> rAddm = SHE.addList(r, encX, paramters);
				
				rList.put(key, r);
				rAddmList.put(key, rAddm);
			}
			
		}
		
		HashMap<Integer, ArrayList<BigInteger>> result = KDTree.queryResultRecovery(rList, rAddmList, paramters);
		
		
		return result;
	}

}
