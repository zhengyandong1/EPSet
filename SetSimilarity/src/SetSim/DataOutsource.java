package SetSim;
import java.math.BigInteger;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Set;

import SetSim.KDTree.EncNode;
import SetSim.KDTree.Node;


public class DataOutsource {
	
	public static int k = MainFinal.k;
	public static int upBoundLeafNode = MainFinal.upBoundLeafNode;
	public static BigInteger cipherZero1 = MainFinal.cipherZero1;
	public static BigInteger cipherZero2 = MainFinal.cipherZero2; 
	public static BigInteger cipherMinusOne = MainFinal.cipherMinusOne; 
	public static SHE paramters = MainFinal.paramters;
	
	public static ArrayList<Integer> vectorToSet(double[] vector) {
		
		
		
		ArrayList<Integer> set = new ArrayList<Integer>();
		
		for(int i = 0; i < vector.length; i++) {
			
			if(vector[i] == 1) {
				set.add(i);
			}
		}
			
		
		return set;
		
	}

	
	public static HashMap<Integer, ArrayList<Integer>> vectorMapToSetList(HashMap<Integer, double[]> dataList) {
		
		
		HashMap<Integer, ArrayList<Integer>> setList = new HashMap<Integer, ArrayList<Integer>>();
		
		for (Iterator iterator = dataList.keySet().iterator(); iterator.hasNext();) {
			
			int key = (int) iterator.next();
			double[] vector = dataList.get(key);
			
			ArrayList<Integer> set = new ArrayList<Integer>();
			
			for(int i = 0; i < vector.length; i++) {
				
				if(vector[i] == 1) {
					set.add(i);
				}
			}
			
			setList.put(key, set);
			
		}
		
		return setList;
		
	}

	
	public static HashMap<Integer, ArrayList<Integer>> vectorListToSetList(ArrayList<double[]> pivotList) {
		
		
		HashMap<Integer, ArrayList<Integer>> setList = new HashMap<Integer, ArrayList<Integer>>();
		
		for(int i = 0; i < pivotList.size(); i++) {
			
			double[] vector = pivotList.get(i);
			
			ArrayList<Integer> pivot = new ArrayList<Integer>();
			
			for(int j = 0; j < vector.length; j++) {
				if(vector[j] == 1) {
					pivot.add(j);
				}
			}
			
			setList.put(i, pivot);
			
		}
		
		return setList;
		
		
	}
	
	public static EncNode dataOutsourcing(HashMap<Integer, double[]> dataList, ArrayList<double[]> pivotList) {
		
		int n = dataList.size();
		
		// choose pivots
		Set<Integer> pivotLocs = choosePivots(dataList, k);

		for (Iterator iterator = pivotLocs.iterator(); iterator.hasNext();) {
			
			int key = (int) iterator.next();
			pivotList.add(dataList.get(key));
			
		}
		
		HashMap<Integer, double[]> distanceList = new HashMap<Integer, double[]>();
		
		for(int i = 0; i < dataList.size(); i++) {

			double[] distances = DistanceComPivot(dataList.get(i), pivotList);
			distanceList.put(i, distances);
			
		}
		
		Node root = KDTree.TreeBuild(distanceList, upBoundLeafNode, k);
		
		HashMap<Integer, ArrayList<Integer>> setList = vectorMapToSetList(dataList);
		
		EncNode encryptRoot = KDTree.encryptTree(setList, root, paramters);
		
		return encryptRoot;
		
	}
	
	
	public static Node treeBuildPlain(HashMap<Integer, double[]> dataList, ArrayList<double[]> pivotList) {
		
		int n = dataList.size();
		
		// choose pivots
		Set<Integer> pivotLocs = choosePivots(dataList, k);

		for (Iterator iterator = pivotLocs.iterator(); iterator.hasNext();) {
			
			int key = (int) iterator.next();
			pivotList.add(dataList.get(key));
			
		}
		
		HashMap<Integer, double[]> distanceList = new HashMap<Integer, double[]>();
		
		for(int i = 0; i < dataList.size(); i++) {

			double[] distances = DistanceComPivot(dataList.get(i), pivotList);
			distanceList.put(i, distances);
			
		}
		
		Node root = KDTree.TreeBuild(distanceList, upBoundLeafNode, k);
		
		return root;
	}
	
	public static Set<Integer> choosePivots(HashMap<Integer, double[]> dataLists, int k) {
		
		Set<Integer> pivotLocs = new HashSet<Integer>();
		int n = dataLists.size();
		
		// one pivot
		
		pivotLocs.add(0);
		
		for(int i = 0; i < k - 1; i++) {
			
			double sim = 0.0;
			int loc = 0;
			for(int j = 0; j < n; j++) {
				if(!pivotLocs.contains(j)) {
					sim = Sim(dataLists.get(j), dataLists, pivotLocs);
					loc = j;
					break;
				}
			}
			
			for(int j = loc + 1; j < n; j++) {
				if(!pivotLocs.contains(j)) {
					double tempSim = Sim(dataLists.get(j), dataLists, pivotLocs);
					if(tempSim < sim) {
						sim = tempSim;
						loc = j;
					}
				}
			}
			
			pivotLocs.add(loc);
			
		}
		
		return pivotLocs;
	}
	
	public static double Sim(double[] a, HashMap<Integer, double[]> dataLists, Set<Integer> pivotLocs) {
		
		double sim = 0.0;
		
		for (Iterator iterator = pivotLocs.iterator(); iterator.hasNext();) {
			Integer loc = (Integer) iterator.next();
			sim = sim + Jaccard(a, dataLists.get(loc));
		}
		
		return sim;
	}
	
	public static double Jaccard(double[] a, double[] b) {
		
		int n = a.length;
		int intersection = 0;
		int union = 0;
		
		for(int i = 0; i < n; i++) {
			if(a[i] == 1 && b[i] == 1) {
				intersection = intersection + 1;
				union = union + 1;
			}else if(a[i] == 1 || b[i] == 1){
				union = union + 1;
			}
		}
		return intersection*1.0/union;
	}
	
	
	public static double[] DistanceComPivot(double[] data, ArrayList<double[]> pivotList) {
		
		double[] distances = new double[pivotList.size()];
		
		for(int i = 0; i < distances.length; i++) {
			distances[i] = 1 - Jaccard(data, pivotList.get(i));
		}
		
		return distances;
	}
	
	

}
