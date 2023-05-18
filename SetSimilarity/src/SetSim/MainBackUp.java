package SetSim;

import java.math.BigInteger;
import java.util.ArrayList;
import java.util.HashMap;

import SetSim.KDTree.EncNode;
import SetSim.Trapdoor.FilterTrapdoor;
import SetSim.Trapdoor.RefinementTrapdoor;

public class MainBackUp {
	
		
	// the number of pivots
	public static int k = 4;
	public static int scale = 1000;
	public static int upBoundLeafNode = 5; 
	
	// security parameters
	public static int k0 = 2048;
	public static int k1 = 40;
	public static int k2 = 160;
	
	public static SHE paramters;
	public static BigInteger cipherZero1;
	public static BigInteger cipherZero2; 
	public static BigInteger cipherMinusOne;
	
	public static void keyGen(int k0, int k1, int k2) {
		
		paramters = SHE.keyGen(k0, k1, k2);
		
		cipherZero1 = SHE.enc(0, paramters);
		cipherZero2 = SHE.enc(0, paramters);
		cipherMinusOne = SHE.enc(-1, paramters);
		
	}
	
	public static double dataOutsourcing(int n, int l) throws Exception {
		
		int cycle = 50;
		double t = 0.0;
		
		HashMap<Integer, double[]> dataList = Data.ReadDataJeterWithoutDup(n, l);
		ArrayList<double[]> pivotList = new ArrayList<double[]>();
		
		for(int v = 0; v < cycle; v++) {
			
			double t1 = System.currentTimeMillis();
			EncNode encNode = DataOutsource.dataOutsourcing(dataList, pivotList);
			double t2 = System.currentTimeMillis();
			
			if(v >= cycle*0.8) {
				t = t + (t2 - t1);
			}
		}
		
		return t*5/cycle;
		
	}
	
	
	
	
	
	public static double[] queryProcessing(int n, int l) throws Exception {
		
		int cycle = 50;
		int tau1 = 90;
		int tau2 = 100;
		
		double[] t = new double[2];
		
		HashMap<Integer, double[]> dataList = Data.ReadDataJeterWithoutDup(n, l);
		
		ArrayList<double[]> pivotList = new ArrayList<double[]>();
		EncNode encNode = DataOutsource.dataOutsourcing(dataList, pivotList);
		
		HashMap<Integer, ArrayList<Integer>> pivotHashMap = DataOutsource.vectorListToSetList(pivotList);
		
		for(int v = 0; v < cycle; v++) {
			
			int index = (int) (n*Math.random());
			double[] qvector = dataList.get(index);
			
			ArrayList<Integer> querySet = DataOutsource.vectorToSet(qvector);
			
			FilterTrapdoor filterTrapdoor = Trapdoor.genFilterTrapdoor(querySet, pivotHashMap, tau1, tau2, paramters);
			RefinementTrapdoor refinementTrapdoor = Trapdoor.genRefinementTrapdoor(querySet, tau1, tau2, paramters);
			
			HashMap<Integer, ArrayList<BigInteger>> rList = new HashMap<Integer, ArrayList<BigInteger>>();
			HashMap<Integer, ArrayList<BigInteger>> rAddmList = new HashMap<Integer, ArrayList<BigInteger>>();
			
			double t1 = System.currentTimeMillis();
			KDTree.SearchTree(encNode, filterTrapdoor, refinementTrapdoor, paramters, rList, rAddmList);
			double t2 = System.currentTimeMillis();
			
//			HashMap<Integer, ArrayList<BigInteger>> result = KDTree.queryResultRecovery(rList, rAddmList, paramters);
//			HashMap<Integer, ArrayList<Integer>> resultInteger = integerListToBigIntegerList(result);
//			
//			System.out.println("query results of our scheme");
//			
//			for (Integer key : result.keySet()) {
//				
//				print(result.get(key));
//				System.out.println(Trapdoor.jaccardComp(querySet, resultInteger.get(key)));
//				
//			}
			
			// query by Naive
//			System.out.println("query results by Naive");
			HashMap<Integer, ArrayList<BigInteger>> rNaiveList = new HashMap<Integer, ArrayList<BigInteger>>();
			HashMap<Integer, ArrayList<BigInteger>> rNaiveAddmList = new HashMap<Integer, ArrayList<BigInteger>>();
			
			HashMap<Integer, ArrayList<BigInteger>> encDataNaive = Naive.encDataSetNaive(n, l);
			
			double t3 = System.currentTimeMillis();
			HashMap<Integer, ArrayList<BigInteger>> result = Naive.queryNaive(encDataNaive, refinementTrapdoor, rNaiveList, rNaiveAddmList);
			double t4 = System.currentTimeMillis();
			
			if(v >= cycle*0.8) {
				t[0] = t[0] + (t2 - t1);
				t[1] = t[1] + (t3 - t4);
			}
			
			
			System.out.printf("our scheme: %f; Naive scheme: %f", (t2 - t1), (t4 - t3));
			
//			resultInteger = integerListToBigIntegerList(result);
//			
//			
//			for (Integer key : result.keySet()) {
//				
//				print(result.get(key));
//				System.out.println(Trapdoor.jaccardComp(querySet, resultInteger.get(key)));
//				
//			}
//			
//			System.out.println("query results over plaintexts");
//			HashMap<Integer, ArrayList<Integer>> resultPlain = QueryPlain.queryPlain(dataList, querySet, tau1, tau2);
//			
//			
//			for (Integer key : resultPlain.keySet()) {
//				
//				printInteger(resultPlain.get(key));
//				System.out.println(Trapdoor.jaccardComp(querySet, resultPlain.get(key)));
//				
//			}
				
		}
		
		System.out.printf("our scheme: %f; Naive scheme: %f", (t[0]*5/cycle), t[1]*5/cycle);
		
		return t;
		
	}
	
	
	public static void dataOutsource(int n, int l) throws Exception {
		
		HashMap<Integer, double[]> dataList = Data.ReadDataJeterWithoutDup(n, l);
		ArrayList<double[]> pivotList = new ArrayList<double[]>();
		
		EncNode encNode = DataOutsource.dataOutsourcing(dataList, pivotList);
		
		HashMap<Integer, ArrayList<Integer>> pivotHashMap = DataOutsource.vectorListToSetList(pivotList);
		
		for(int cycle = 0; cycle <= 10; cycle++) {
			
			int index = (int) (n*Math.random());
			double[] qvector = dataList.get(index);
			
			ArrayList<Integer> querySet = DataOutsource.vectorToSet(qvector);
			
			System.out.println("query set with index:" + index);
//			printInteger(querySet);
			
			int tau1 = 80;
			int tau2 = 100;
			
			System.out.printf("tau1 = %d, tau2 = %d", tau1, tau2);
			
			FilterTrapdoor filterTrapdoor = Trapdoor.genFilterTrapdoor(querySet, pivotHashMap, tau1, tau2, paramters);
			RefinementTrapdoor refinementTrapdoor = Trapdoor.genRefinementTrapdoor(querySet, tau1, tau2, paramters);
			
			HashMap<Integer, ArrayList<BigInteger>> rList = new HashMap<Integer, ArrayList<BigInteger>>();
			HashMap<Integer, ArrayList<BigInteger>> rAddmList = new HashMap<Integer, ArrayList<BigInteger>>();
			
			
			double t1 = System.currentTimeMillis();
			KDTree.SearchTree(encNode, filterTrapdoor, refinementTrapdoor, paramters, rList, rAddmList);
			double t2 = System.currentTimeMillis();
			
			System.out.println("query time:" + (t2 - t1));
			
			HashMap<Integer, ArrayList<BigInteger>> result = KDTree.queryResultRecovery(rList, rAddmList, paramters);
			HashMap<Integer, ArrayList<Integer>> resultInteger = integerListToBigIntegerList(result);
			
			System.out.println("query results of our scheme");
			
			for (Integer key : result.keySet()) {
				
				print(result.get(key));
				System.out.println(Trapdoor.jaccardComp(querySet, resultInteger.get(key)));
				
			}
			
			// query by Naive
			System.out.println("query results by Naive");
			HashMap<Integer, ArrayList<BigInteger>> rNaiveList = new HashMap<Integer, ArrayList<BigInteger>>();
			HashMap<Integer, ArrayList<BigInteger>> rNaiveAddmList = new HashMap<Integer, ArrayList<BigInteger>>();
			
			HashMap<Integer, ArrayList<BigInteger>> encDataNaive = Naive.encDataSetNaive(n, l);
			result = Naive.queryNaive(encDataNaive, refinementTrapdoor, rNaiveList, rNaiveAddmList);
			
			resultInteger = integerListToBigIntegerList(result);
			
			
			for (Integer key : result.keySet()) {
				
				print(result.get(key));
				System.out.println(Trapdoor.jaccardComp(querySet, resultInteger.get(key)));
				
			}
			
			
//			System.out.println("query results over plaintexts");
//			HashMap<Integer, ArrayList<Integer>> resultPlain = QueryPlain.queryPlain(dataList, querySet, tau1, tau2);
//			
//			for (Integer key : resultPlain.keySet()) {
//				
//				printInteger(resultPlain.get(key));
//				System.out.println(Trapdoor.jaccardComp(querySet, resultPlain.get(key)));
//				
//			}
			
			
//			System.out.println("query result over plaintext tree");
//			
//			double tau = tau1*1.0/tau2;
//			
//			double[] leftDisVector = new double[k];
//			double[] rightDisVector = new double[k];
//			
//			
//			for (int i = 0; i < pivotHashMap.keySet().size(); i++) {
//				
//				
//				double sim = Trapdoor.jaccardComp(querySet, pivotHashMap.get(i));
//				
//				double dis = 1 - sim;
//				
//				leftDisVector[i] = (int)((dis - (1 - tau))*scale);
//				rightDisVector[i] = (int) ((dis + (1 - tau))*scale);
//				
//			}
//			
//			
//			Node root = DataOutsource.treeBuildPlain(dataList, pivotList);
//			
//			KDTree.searchPlain(root, querySet, leftDisVector, rightDisVector, tau1, tau2);
//		
		}
		
	}
	
	
	public static void print(ArrayList<BigInteger> list) {
		for(int i = 0; i < list.size(); i++) {
			System.out.printf("%d,", list.get(i));
		}
		System.out.println();
	}
	
	public static void printInteger(ArrayList<Integer> list) {
		for(int i = 0; i < list.size(); i++) {
			System.out.printf("%d,", list.get(i));
		}
		System.out.println();
	}
	
	
	public static HashMap<Integer, ArrayList<Integer>> integerListToBigIntegerList(HashMap<Integer, ArrayList<BigInteger>> result) {
		
		HashMap<Integer, ArrayList<Integer>> resultIntegerList = new HashMap<Integer, ArrayList<Integer>>();
		
		for (Integer key : result.keySet()) {
			
			ArrayList<Integer> integerList = new ArrayList<Integer>();
			
			for(int i = 0; i < result.get(key).size(); i++) {
				integerList.add(result.get(key).get(i).intValue());
			}
			
			resultIntegerList.put(key, integerList);
			
		}
		
		return resultIntegerList;
		
	}
	
	public static void main(String[] args) throws Exception {
		
		int n = 10000;
		int l = 50;
		
		keyGen(k0, k1, k2);
		
		System.out.println("start MainFinal");
		
//		dataOutsource(n, l);
		queryProcessing(n, l);
		
		
		System.out.println("end");
		
	}


}
