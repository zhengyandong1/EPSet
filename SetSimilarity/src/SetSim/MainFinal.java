package SetSim;

import java.io.FileOutputStream;
import java.io.PrintStream;
import java.math.BigInteger;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Random;

import SetSim.KDTree.EncNode;
import SetSim.KDTree.InterNode;
import SetSim.KDTree.Node;
import SetSim.Le.Ciphertext;
import SetSim.Le.TrapdoorLe;
import SetSim.Trapdoor.FilterTrapdoor;
import SetSim.Trapdoor.RefinementTrapdoor;

public class MainFinal {
	
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
	
	// Computational costs of internal nodes encryption
	public static double internalNodeEnc(int k) {
		
		double t = 0;
		int cycle = 1000;
		Random rand = new Random();
		
		for(int j = 0; j < cycle; j++) {
			
			double t1 =  System.nanoTime();
			int c = rand.nextInt(k);
			int val = rand.nextInt();
			
			BigInteger[] b_c = new BigInteger[k];
			
			for(int i = 0; i < b_c.length; i++) {
				if(i== c) {
					b_c[i] = SHE.enc(1, paramters);
				}else {
					b_c[i] = SHE.enc(0, paramters);
				}
			}
			BigInteger minusVal = SHE.enc(-val, paramters);
			double t2 = System.nanoTime();
			
			if(j >= cycle*0.8) {
				t = t + (t2 - t1);
			}
		}
		
		return t*5/cycle/1000000;
	}
	
	public static double leafNodeEnc(int sizeSet) {
		
		
		double t = 0;
		int cycle = 1000;
		Random rand = new Random();
		
		ArrayList<Integer> record = new ArrayList<>();
		
		for(int i = 0; i < sizeSet; i++) {
			record.add(rand.nextInt());
		}
		
		for(int j = 0; j < cycle; j++) {
			double t1 = System.nanoTime();
			ArrayList<BigInteger> encRecord = SHE.encIntegerList(record, paramters);
			double t2 = System.nanoTime();
			
			if(j >= cycle*0.8) {
				t = t + (t2 - t1);
			}
		}
		
		t = t*5/cycle/1000000;
		
		return t;
		
	}
	
	// Computational costs of data outsourcing
	public static double dataOutsourcing(int n, int l, int k) throws Exception {
		
		MainFinal.k = k;
		
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
	
	// Computational costs of trapdoor generation
	public static double[] trapdoorgenearion(int k, int sizeSet) {
		
		MainFinal.k = k;
		double[] t = new double[2];
		int cycle = 10000;
		
		Random rand = new Random();
		ArrayList<Integer> qSet = new ArrayList<Integer>();
		for(int i = 0; i < sizeSet; i++) {
			qSet.add(rand.nextInt());
		}
		
		HashMap<Integer, ArrayList<Integer>> pivots = new HashMap<Integer, ArrayList<Integer>>();
		
		for(int i = 0; i < k; i++) {
			pivots.put(i, qSet);
		}
		
		int tau1 = 90;
		int tau2 = 100;
	
		for(int i = 0; i < cycle; i++) {
			
			// generate filter trapdoor
			double t1 = System.nanoTime();
			FilterTrapdoor filterTrapdoor = Trapdoor.genFilterTrapdoor(qSet, pivots, tau1, tau2, paramters);
			double t2 = System.nanoTime();
			
			// generate refinement trapdoor
			double t3 = System.nanoTime();
			RefinementTrapdoor refinementTrapdoor = Trapdoor.genRefinementTrapdoor(qSet, tau1, tau2, paramters);
			double t4 = System.nanoTime();
			
			if(i >= cycle*0.8) {
				t[0] = t[0] + (t2 - t1);
				t[1] = t[1] + (t4 - t3);
			}
			
		}
		
		t[0] = t[0]*5/cycle/1000000;
		t[1] = t[1]*5/cycle/1000000;
		
		return t;
		
	}
	
	public static double leafnodeQuery(int sizeSet1, int sizeSet2) {
		
		double t = 0;
		int cycle = 100;
		Random rand = new Random();
		
		ArrayList<Integer> record = new ArrayList<>();
		for(int i = 0; i < sizeSet1; i++) {
			record.add(rand.nextInt());
		}
		ArrayList<BigInteger> encRecord = SHE.encIntegerList(record, paramters);
		
		int tau1 = 90;
		int tau2 = 100;
		
		ArrayList<Integer> qSet = new ArrayList<>();
		for(int i = 0; i < sizeSet2; i++) {
			qSet.add(rand.nextInt());
		}
		
		RefinementTrapdoor refinementTrapdoor = Trapdoor.genRefinementTrapdoor(qSet, tau1, tau2, paramters);
		
		for(int j = 0; j < cycle; j++) {
			double t1 = System.nanoTime();
			boolean result = KDTree.leafEval(encRecord, refinementTrapdoor, paramters);
			double t2 = System.nanoTime();
			
			if(j >= cycle*0.8) {
				t = t + (t2 - t1);
			}
		}
		
		t = t*5/cycle/1000000;
		
		return t;
	}
	
	// computational costs of query processing
	public static double queryProcessing(int n, int l, int k, int tau1, int tau2) throws Exception {
		
		
		int cycle = 20;
		double t = 0;
		
		HashMap<Integer, double[]> dataList = Data.ReadDataJeterWithoutDup(n, l);
		ArrayList<double[]> pivotList = new ArrayList<double[]>();
		
		// data encryption using our scheme
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
			
			// query processing by our scheme 
			double t1 = System.currentTimeMillis();
			KDTree.SearchTree(encNode, filterTrapdoor, refinementTrapdoor, paramters, rList, rAddmList);
			double t2 = System.currentTimeMillis();
					
			
			t = t + (t2 - t1);
				
		}
		
		return t/cycle;
		
	}
	
	public static double[] queryProcessingComp(int n, int l, int k, int tau1, int tau2, int overElementSize) throws Exception {
		
		int cycle = 20;
		double[] t = new double[2];
		
		HashMap<Integer, double[]> dataList = Data.ReadDataJeterWithoutDup(n, l);
		ArrayList<double[]> pivotList = new ArrayList<double[]>();
		
		// data encryption using our scheme
		EncNode encNode = DataOutsource.dataOutsourcing(dataList, pivotList);
		HashMap<Integer, ArrayList<Integer>> pivotHashMap = DataOutsource.vectorListToSetList(pivotList);
		
		// data encryption naive
		HashMap<Integer, ArrayList<BigInteger>> encDataNaive = Naive.encDataSetNaive(n, l);
		
		
		for(int v = 0; v < cycle; v++) {
			
			int index = (int) (n*Math.random());
			double[] qvector = dataList.get(index);
			
			ArrayList<Integer> querySet = DataOutsource.vectorToSet(qvector);
			
			FilterTrapdoor filterTrapdoor = Trapdoor.genFilterTrapdoor(querySet, pivotHashMap, tau1, tau2, paramters);
			RefinementTrapdoor refinementTrapdoor = Trapdoor.genRefinementTrapdoor(querySet, tau1, tau2, paramters);
			
			HashMap<Integer, ArrayList<BigInteger>> rList = new HashMap<Integer, ArrayList<BigInteger>>();
			HashMap<Integer, ArrayList<BigInteger>> rAddmList = new HashMap<Integer, ArrayList<BigInteger>>();
			
			// query processing by our scheme 
			double t1 = System.currentTimeMillis();
			KDTree.SearchTree(encNode, filterTrapdoor, refinementTrapdoor, paramters, rList, rAddmList);
			double t2 = System.currentTimeMillis();
			
			// query processing by naive
			HashMap<Integer, ArrayList<BigInteger>> rNaiveList = new HashMap<Integer, ArrayList<BigInteger>>();
			HashMap<Integer, ArrayList<BigInteger>> rNaiveAddmList = new HashMap<Integer, ArrayList<BigInteger>>();
			
			double t3 = System.currentTimeMillis();
			HashMap<Integer, ArrayList<BigInteger>> result = Naive.queryNaive(encDataNaive, refinementTrapdoor, rNaiveList, rNaiveAddmList);
			double t4 = System.currentTimeMillis();
			
			
			
			if(v >= cycle*0.8) {
				t[0] = t[0] + (t2 - t1);
				t[1] = t[1] + (t4 - t3);
			}
				
		}
		
		t[0] = t[0]*5/cycle;
		t[1] = t[1]*5/cycle;
		
		return t;
		
	}
	
//	public static double queryProcessingNaive(int n, int l, int tau1, int tau2) throws Exception {
//		
//		double t = 0;
//		int cycle = 20;
//		HashMap<Integer, double[]> dataList = Data.ReadDataJeterWithoutDup(n, l);
//		HashMap<Integer, ArrayList<BigInteger>> encDataNaive = Naive.encDataSetNaive(n, l);
//		
//		for(int v = 0; v < cycle; v++) {
//			
//			int index = (int) (n*Math.random());
//			double[] qvector = dataList.get(index);
//			
//			ArrayList<Integer> querySet = DataOutsource.vectorToSet(qvector);
//			
//			RefinementTrapdoor refinementTrapdoor = Trapdoor.genRefinementTrapdoor(querySet, tau1, tau2, paramters);
//			
//			HashMap<Integer, ArrayList<BigInteger>> rNaiveList = new HashMap<Integer, ArrayList<BigInteger>>();
//			HashMap<Integer, ArrayList<BigInteger>> rNaiveAddmList = new HashMap<Integer, ArrayList<BigInteger>>();
//			
//			double t1 = System.currentTimeMillis();
//			HashMap<Integer, ArrayList<BigInteger>> result = Naive.queryNaive(encDataNaive, refinementTrapdoor, rNaiveList, rNaiveAddmList);
//			double t2 = System.currentTimeMillis();
//			
//			if(v >= cycle*0.8) {
//				t = t + (t2 - t1);
//			}
//		}
//		
//		return t*5/cycle;
//		
//	}
//	
//	public static double queryProcessingLe(int n, int overElementSize) throws Exception {
//		
//		int l = 100;
//		int cycle = 20;
//		double t = 0.0;
//		
//		HashMap<Integer, double[]> dataList = Data.ReadDataJeterWithoutDup(n, l);
//		ArrayList<Ciphertext> encDataset = Le.encDataLe(n, overElementSize);
//		
//		TrapdoorLe trapdoorLe = Le.trapdoorGen(overElementSize);
//		
//		for(int v = 0; v < cycle; v++) {
//			
//			double t1 = System.currentTimeMillis();
//			Le.queryLe(encDataset, trapdoorLe);
//			double t2 = System.currentTimeMillis();
//			
//			if(v >= 0.8*cycle) {
//				t = t + (t2 - t1);
//			}
//		}
//		
//		return t*5/cycle;
//		
//	}
	
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
	
	// data outsourcing evaluation
	public static void dataOutsourcingEval() throws Exception {
		
		int[] n = {4000, 8000, 12000, 16000,  20000};
		int[] k = {4, 6, 8, 10};
		int[] setSize = {5, 10, 15, 20, 25};
		
		int l = 100;
		
		// dataset encryption with n and k
		System.out.println("Computational costs of dataset encryption with n and k");
		System.out.printf("\t");
		for(int i = 0; i < n.length; i++) {
			System.out.printf("n = %d\t", n[i]);
		}
		System.out.println();
		
		for(int j = 0; j < k.length; j++) {
			System.out.printf("k = %d\t", k[j]);
			for(int i = 0; i < n.length; i++) {
				System.out.printf("%f\t", dataOutsourcing(n[i], l, k[j]));
			}
			
			System.out.println();
		}
		
		System.out.println();
		
	}
	
	public static void trapdoorGenEval() {
		
		int[] k = {4, 6, 8, 10};
		int[] setSize = {5, 10, 15, 20, 25};
		
		// trapdoor generation with k and setSize
		System.out.println("Computational costs of trapdoor encryption with k and setSize");
		
		System.out.printf("\t");
		for(int i = 0; i < setSize.length; i++) {
			System.out.printf("setSize = %d\t", setSize[i]);
		}
		System.out.println();
		
		for(int j = 0; j < k.length; j++) {
			System.out.printf("k = %d\t", k[j]);
			for(int i = 0; i < setSize.length; i++) {
				double[] t = trapdoorgenearion(k[j], setSize[i]);
				System.out.printf("%f\t", (t[0] + t[1]));
			}
			
			System.out.println();
		}
	}
	
	public static void queryProcessingEval() throws Exception {
		
		int[] n = {4000, 8000, 12000, 16000,  20000};
		int[] k = {4, 6, 8};
		int[] tau = {90, 92, 94, 96, 98};
		int[] setSize = {5, 10, 15, 20, 25};
		
		int nVal = 20000;
		int kVal = 6;
		int tauVal = 95;
		
		int l = 100;
		
		// query processing with n and k
		System.out.println("Computational costs of query processing with n and k");
		System.out.printf("\t");
		for(int i = 0; i < n.length; i++) {
			System.out.printf("n = %d\t", n[i]);
		}
		System.out.println();
		
		for(int j = 0; j < k.length; j++) {
			System.out.printf("k = %d\t", k[j]);
			for(int i = 0; i < n.length; i++) {
				double t = queryProcessing(n[i], l, k[j], tauVal, 100);
				System.out.printf("%f\t", t);
			}
			
			System.out.println();
		}
		
		System.out.println();
		
		// query processing with tau and k
//		System.out.println("Computational costs of query processing with tau and k");
//		System.out.printf("\t");
//		for(int i = 0; i < tau.length; i++) {
//			System.out.printf("tau = %d\t", tau[i]);
//		}
//		System.out.println();
//		
//		for(int j = 0; j < k.length; j++) {
//			System.out.printf("k = %d\t", k[j]);
//			for(int i = 0; i < tau.length; i++) {
//				double t = queryProcessing(nVal, l, k[j], tau[i], 100);
//				System.out.printf("%f\t", t);
//			}
//			System.out.println();
//		}
//		
//		System.out.println();
		
	}
	
	public static void leafNodeQueryEval() {
		
		int overallSetSize = 35;
		
		for(int i = 0; i < overallSetSize; i++) {
			for(int j = 0; j < overallSetSize; j++) {
				System.out.printf("%f\t", leafnodeQuery(i, j));
			}
			System.out.println();
		}
		
	}
	
	public static void main(String[] args) throws Exception {
		
		keyGen(k0, k1, k2);

		FileOutputStream bos = new FileOutputStream("src/output.txt");
		System.setOut(new PrintStream(bos));
		
//		Computational costs of our scheme
		
		System.out.println("Our scheme: data outsourcing with k and tau");
		queryProcessingEval();
		
//		Evaluate Le scheme
//		System.out.println("Le: data outsourcing with n and Uset");
//		Le.encDataMultiEval();
//		System.out.println("Le: trapdoor generation with Uset");
//		Le.trapdoorGenEval();
//		System.out.println("Le: query processing with n and Uset");
//		Le.queryLeEval();
//
//		System.out.println("Our scheme: data outsoucring with n and k");
//		dataOutsourcingEval();
//		System.out.println("Our scheme: trapdoor generation with k and USet");
//		trapdoorGenEval();
		
		
	}
	
	

}
