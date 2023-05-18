package SetSim;

import java.math.BigInteger;
import java.security.SecureRandom;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Stack;


public class KDTree {
	
	public static int k = MainFinal.k;
	public static int scale = MainFinal.scale; 
	
	public static BigInteger cipherZero1 = MainFinal.cipherZero1;
	public static BigInteger cipherZero2 = MainFinal.cipherZero2; 
	public static BigInteger cipherMinusOne = MainFinal.cipherMinusOne; 
	
	public static SecureRandom rnd = new SecureRandom();
	
	public class Node{
		public boolean isLeafNode;
	}
	
	public class InterNode extends Node{
		
		public int c;
		public double val;
		public Node left;
		public Node right;
		
		public InterNode(int c, double val, Node left, Node right) {
			
			this.c = c;
			this.val = val;
			this.left = left;
			this.right = right;
			this.isLeafNode = false;
			
		}
	}
	
	public class LeafNode extends Node {

		public HashMap<Integer, double[]> setRecords;

		public LeafNode(HashMap<Integer, double[]> setRecords) {
			
			this.setRecords = setRecords;
			this.isLeafNode = true;
			
		}
		
	}
	
	public class EncNode{
		public boolean isLeafNode;
	}
	
	public class EncInternal extends EncNode{
		
		public BigInteger[] b_c;
		public BigInteger minusVal;
		public EncNode left;
		public EncNode right;
		
		public EncInternal(BigInteger[] b_c, BigInteger minusVal, KDTree.EncNode left, KDTree.EncNode right) {
			
			this.b_c = b_c;
			this.minusVal = minusVal;
			this.left = left;
			this.right = right;
			this.isLeafNode = false;
			
		}
		
	}
	
	public class EncLeaf extends EncNode{
		
		public HashMap<Integer, ArrayList<BigInteger>> encSetRecords;

		public EncLeaf(HashMap<Integer, ArrayList<BigInteger>> encSetRecords) {
			
			this.encSetRecords = encSetRecords;
			this.isLeafNode = true;
			
		}
		
	}

	
	public static KDTree kdTree = new KDTree();
	
	public static Node TreeBuild(HashMap<Integer, double[]> dataList, int upBoundLeafNode, int k) {
		
		int N = dataList.size();
		
		if(N <= upBoundLeafNode) {
			
			LeafNode leafNode = kdTree.new LeafNode(dataList);
			return leafNode;
		}
		
		int c = chooseCD(dataList, k);
		
		double[] dataInC = new double[N];
		
		int i = 0;
		for (Iterator iterator = dataList.keySet().iterator(); iterator.hasNext();) {
			int key = (int) iterator.next();
			dataInC[i] = dataList.get(key)[c];
			i++;
		}
		
		Arrays.sort(dataInC);
		double max = dataInC[dataInC.length - 1];
		double val = Median(dataInC);
		
		if(val == max) {
			LeafNode leafNode = kdTree.new LeafNode(dataList);
			return leafNode;
		}
		
		HashMap<Integer, double[]> leftDataList = new HashMap<Integer, double[]>();
		HashMap<Integer, double[]> rightDataList = new HashMap<Integer, double[]>();
		
		
		HashMap<Integer, ArrayList<Integer>> leftSetRecords = new HashMap<Integer, ArrayList<Integer>>();
		HashMap<Integer, ArrayList<Integer>> rightSetRecords = new HashMap<Integer, ArrayList<Integer>>();
		
		for (Iterator iterator = dataList.keySet().iterator(); iterator.hasNext();) {
			
			int key = (int) iterator.next();
			if(dataList.get(key)[c] <= val) {
				leftDataList.put(key, dataList.get(key));
			}else {
				rightDataList.put(key, dataList.get(key));
			}
			
		} 
		
		Node left = TreeBuild(leftDataList, upBoundLeafNode, k);
		Node right = TreeBuild(rightDataList, upBoundLeafNode, k);
		
		InterNode interNode = kdTree.new InterNode(c, val, left, right);
		
		return interNode;
		
	}

	public static void searchPlain(Node root, ArrayList<Integer> qSet, double[] leftRange, double[] rightRange, int tau1, int tau2) {
		
		HashMap<Integer, ArrayList<Integer>> result = new HashMap<Integer, ArrayList<Integer>>();
		
		Stack<Node> S = new Stack<Node>();
		S.push(root);
		
		while(S.size() > 0) {
			
			Node node = S.pop();
			
			if(node.isLeafNode == true) {
				
				LeafNode leafNode = (LeafNode) node;
				
				HashMap<Integer, double[]> set = leafNode.setRecords;
				HashMap<Integer, ArrayList<Integer>> setInteger = DataOutsource.vectorMapToSetList(set);
				
				for (Iterator iterator = setInteger.keySet().iterator(); iterator.hasNext();) {
					
					int key = (int) iterator.next();
					if(Trapdoor.jaccardComp(qSet, setInteger.get(key)) >= tau1*1.0/tau2) {
						MainFinal.printInteger(setInteger.get(key));
					}
					
				}
			}else{
				
				InterNode interNode = (InterNode) node;
				
				int c = interNode.c;
				double val = interNode.val;
				
				System.out.println(c);
				
				if(val <= leftRange[c]) {
					S.push(interNode.right);
				}else if (val > rightRange[c]) {
					S.push(interNode.left);
				}else {
					S.push(interNode.right);
					S.push(interNode.left);
				}
				
			}
		}
		
	}
	
	public static EncNode encryptTree(HashMap<Integer, ArrayList<Integer>> setRecords, Node node, SHE paramters) {
		
		if(node.isLeafNode == true) {
			
			LeafNode leafNode = (LeafNode) node;
			
			HashMap<Integer, ArrayList<BigInteger>> encSetRecords = new HashMap<Integer, ArrayList<BigInteger>>();
						
			for (Iterator iterator = leafNode.setRecords.keySet().iterator(); iterator.hasNext();) {
				
				int key = (int) iterator.next();
				ArrayList<Integer> record = setRecords.get(key);
				
				ArrayList<BigInteger> encRecord = SHE.encIntegerList(record, paramters);
				encSetRecords.put(key, encRecord);
				
			}
			
			EncLeaf encLeaf = kdTree.new EncLeaf(encSetRecords);
			
			return encLeaf;
		}else{
			
			InterNode interNode = (InterNode) node;
			int c = interNode.c;
			int val = (int) (interNode.val*scale);
			
			BigInteger[] b_c = new BigInteger[k];
			
			for(int i = 0; i < b_c.length; i++) {
				if(i== c) {
					b_c[i] = SHE.enc(1, paramters);
				}else {
					b_c[i] = SHE.enc(0, paramters);
				}
			}
			
			BigInteger minusVal = SHE.enc(-val, paramters);
			
			EncNode left = encryptTree(setRecords, interNode.left, paramters);
			EncNode right = encryptTree(setRecords, interNode.right, paramters);
			
			EncInternal root = kdTree.new EncInternal(b_c, minusVal, left, right);
			
			return root;
			
		}
	} 
	
	public static boolean leafEval(ArrayList<BigInteger> encX, Trapdoor.RefinementTrapdoor refinementTrapdoor, SHE paramters) {
		
		
		BigInteger N = paramters.sheN;
		
		ArrayList<BigInteger> encQ = refinementTrapdoor.encQ;
		BigInteger encTau1 = refinementTrapdoor.encTau1;
		BigInteger encTau2 = refinementTrapdoor.encTau2;
		
		// Computation by S1
		BigInteger[][] matrixXQ = new BigInteger[encQ.size()][encX.size()];
		BigInteger[][] matrixXQply = new BigInteger[encQ.size()][encX.size()];
		
		for(int i = 0; i < encQ.size(); i++) {
			for(int j = 0; j < encX.size(); j++) {
				
				matrixXQ[i][j] = encQ.get(i).multiply(cipherMinusOne).add(encX.get(j)).mod(N);
				matrixXQ[i][j] = SHE.addRandom(matrixXQ[i][j], paramters); 
				
				matrixXQply[i][j] = encX.get(j).multiply(cipherMinusOne).add(encQ.get(i)).mod(N);
				matrixXQply[i][j] = SHE.addRandom(matrixXQply[i][j], paramters); 
			}
		}
		
		// Computation by S2
		BigInteger[][] matrixW = new BigInteger[encQ.size()][encX.size()];
		BigInteger[][] matrixWply = new BigInteger[encQ.size()][encX.size()];
		
		for(int i = 0; i < encQ.size(); i++) {
			for(int j = 0; j < encX.size(); j++) {
				
				matrixXQ[i][j] = SHE.dec(matrixXQ[i][j], paramters);
				matrixXQply[i][j] = SHE.dec(matrixXQply[i][j], paramters);
				
				if(matrixXQ[i][j].compareTo(BigInteger.ZERO) == 1) {
					matrixW[i][j] = SHE.enc(1, paramters);
				}else {
					matrixW[i][j] = SHE.enc(-1, paramters);
				}
				
				if(matrixXQply[i][j].compareTo(BigInteger.ZERO) == 1) {
					matrixWply[i][j] = SHE.enc(1, paramters);
				}else {
					matrixWply[i][j] = SHE.enc(-1, paramters);
				}
			}
		}
		
		// Computation by S1
		BigInteger encSum = BigInteger.ZERO;
		
		for(int i = 0; i < encQ.size(); i++) {
			for(int j = 0; j < encX.size(); j++) {
				encSum = encSum.add(matrixW[i][j]).add(matrixWply[i][j]).mod(N);
			}
		}
		
		BigInteger c1 = encSum.multiply(encTau2).mod(N);
		
		BigInteger c2 = cipherMinusOne.multiply(encTau1).mod(N);
		
		BigInteger c3 =	BigInteger.valueOf(2*encX.size() + 2*encQ.size());
		c3 = c3.add(cipherMinusOne.multiply(encSum)).mod(N);
		
		BigInteger c = c1.add(c2.multiply(c3)).mod(N);
		
		c = SHE.addRandom(c, paramters);
		
		// Computation by S2
		BigInteger m = SHE.dec(c, paramters);
		
		if(m.compareTo(BigInteger.ZERO) == 1) {
			return true;
		}else {
			return false;
		}
		
	}
	
	public static boolean internalLeftEval(BigInteger[] b_c, BigInteger minusVal, Trapdoor.FilterTrapdoor filterTrapdoor, SHE paramters) {
		
		BigInteger N = paramters.sheN;
		
		BigInteger[] leftRange = filterTrapdoor.leftRange;
		
		// Computation at S1
		BigInteger zl = BigInteger.ZERO;
		
		for(int i = 0; i < k; i++) {
			BigInteger cl = leftRange[i].add(minusVal).mod(N);
			cl = SHE.addRandom(cl, paramters);
			cl = b_c[i].multiply(cl).mod(N);
			zl = zl.add(cl).mod(N);
		}
		
		// Computation at S2
		BigInteger zlPlain = SHE.dec(zl, paramters);
		
		if(zlPlain.compareTo(BigInteger.ZERO) == 1) {
			return true;
		}else {
			return false;
		}
		
	}

	
	public static boolean internalRightEval(BigInteger[] b_c, BigInteger minusVal, Trapdoor.FilterTrapdoor filterTrapdoor, SHE paramters) {
		
		BigInteger N = paramters.sheN;
		BigInteger[] rightRange = filterTrapdoor.rightRange;
		
		// Computation at S1
		BigInteger zr = BigInteger.ZERO;
		
		for(int i = 0; i < k; i++) {
			BigInteger cr = rightRange[i].add(minusVal).mod(N);
			cr = SHE.addRandom(cr, paramters);
			cr = b_c[i].multiply(cr).mod(N);
			zr = zr.add(cr).mod(N);
		}
		
		// Computation at S2
		BigInteger zrPlain = SHE.dec(zr, paramters);
		
		if(zrPlain.compareTo(BigInteger.ZERO) == 1) {
			return false;
		}else {
			return true;
		}
		
	}
	
	public static void SearchTree(EncNode root, Trapdoor.FilterTrapdoor filterTrapdoor, Trapdoor.RefinementTrapdoor refinementTrapdoor, SHE paramters, 
			HashMap<Integer, ArrayList<BigInteger>> rList, HashMap<Integer, ArrayList<BigInteger>> rAddmList) {
		
		HashMap<Integer, ArrayList<BigInteger>> result = new HashMap<Integer, ArrayList<BigInteger>>();
		
		Stack<EncNode> S = new Stack<EncNode>();
		S.push(root);
		
		while(S.size() > 0) {
			
			EncNode node = S.pop();
			
			if(node.isLeafNode == true) {
				
				EncLeaf encryptLeafNode = (EncLeaf) node;
				
				HashMap<Integer, ArrayList<BigInteger>> encSetRecords = encryptLeafNode.encSetRecords;
				
				for (Iterator iterator = encSetRecords.keySet().iterator(); iterator.hasNext();) {
					int key = (int) iterator.next();
					if(leafEval(encSetRecords.get(key), refinementTrapdoor, paramters)) {
						
						ArrayList<BigInteger> r = genRandomVector(encSetRecords.get(key).size(), paramters);
						ArrayList<BigInteger> rAddm = SHE.addList(r, encSetRecords.get(key), paramters);
						
						rList.put(key, r);
						rAddmList.put(key, rAddm);
					}
					
				}
			}else{
				
				EncInternal encryptInterNode = (EncInternal) node;
				
				BigInteger[] b_c = encryptInterNode.b_c;
				BigInteger minusVal = encryptInterNode.minusVal;
				
				if(internalLeftEval(b_c, minusVal, filterTrapdoor, paramters)) {
					S.push(encryptInterNode.right);
				}else if (internalRightEval(b_c, minusVal, filterTrapdoor, paramters)) {
					S.push(encryptInterNode.left);
				}else {
					S.push(encryptInterNode.right);
					S.push(encryptInterNode.left);
				}
				
			}
		}
		
	}
	
	public static HashMap<Integer, ArrayList<BigInteger>> queryResultRecovery(HashMap<Integer, ArrayList<BigInteger>> rList, HashMap<Integer, ArrayList<BigInteger>> rAddmList, SHE paramters) {
		
		HashMap<Integer, ArrayList<BigInteger>> resultList = new HashMap<Integer, ArrayList<BigInteger>>(); 
		
		
		for (Iterator iterator = rList.keySet().iterator(); iterator.hasNext();) {
			
			int key = (int) iterator.next();
			ArrayList<BigInteger> r = rList.get(key);
			ArrayList<BigInteger> rAddm = rAddmList.get(key);
			
			ArrayList<BigInteger> result = new ArrayList<BigInteger>();
			
			for(int i = 0; i < r.size(); i++) {
				
				BigInteger rAddmPlain = SHE.dec(rAddm.get(i), paramters);
				result.add(rAddmPlain.subtract(r.get(i)));
			}
			
			resultList.put(key, result);	
		}
		
		return resultList;
		
	}
	
	
	public static ArrayList<BigInteger> genRandomVector(int l, SHE paramters) {
		
		int k1 = paramters.k1;
    	
		ArrayList<BigInteger> randomVector = new ArrayList<BigInteger>();
		
		for(int i = 0; i < l; i++) {
			
			randomVector.add(new BigInteger(k1, rnd));
			
		}
		
		return randomVector;
		
	}
	
	
	public static double Median(double[] dataInC) {	
		
		double median;
		
		int N = dataInC.length;
		if(N%2 == 0) {
			median = 0.5*(dataInC[N/2] + dataInC[N/2 - 1]);
		} else {
			median = dataInC[N/2];
		}
		
		return median;
		
	}
	
	public static int chooseCD(HashMap<Integer, double[]> dataList, int k) {
		
		int c = 0;
		double maxVariance = 0;
		int N = dataList.size();
		
		for(int i = 0; i < k; i++) {
			double[] data = new double[N];
			int j = 0;
			for (Iterator iterator = dataList.keySet().iterator(); iterator.hasNext();) {
				int key = (int) iterator.next();
				data[j] = dataList.get(key)[i];
				j = j + 1;
				
			}
			
			double variance = varianceImperative(data);
			if(variance > maxVariance) {
				c = i;
				maxVariance = variance;
			}
		}
		
		return c;
	}
	
	public static double varianceImperative(double[] population) {
		double average = 0.0;
		for (double p : population) {
			average += p;
		}
		average /= population.length;
 
		double variance = 0.0;
		for (double p : population) {
			variance += (p - average) * (p - average);
		}
		return variance;
	}	
	

}
