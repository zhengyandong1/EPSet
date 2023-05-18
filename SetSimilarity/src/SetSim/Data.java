package SetSim;

import java.io.BufferedReader;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;

public class Data {
	
	public static String filename = "src/jester-data-3.txt";
	public static String splitString = "\t";
	
	
	public static HashMap<Integer, double[]> ReadDataJeterWithoutDup(int n, int l) throws Exception {
		
		HashMap<Integer, double[]> dataList = new HashMap<Integer, double[]>();
		
		String line = null;
	    BufferedReader in = new BufferedReader(new FileReader(filename));
	    
	    int count = 0;
	    while(count < n && (line = in.readLine())!=null)
	    {   
	    	String[] tempStrings = line.split(splitString);
	    	double[] tempDouble = new double[l];
	    	
	    	for(int i = 0; i < l; i++) {
	    		
	    		tempDouble[i] = new Double(tempStrings[i]);
	    		if(tempDouble[i] == 99 || tempDouble[i] <= 0) {
	    			tempDouble[i] = 0;
	    		}else {
	    			tempDouble[i] = 1;
	    		}
	    	}
	    	
			boolean flag = false;
			for (Iterator iterator = dataList.keySet().iterator(); iterator.hasNext();) {
				int key = (int) iterator.next();
				if(Equality(tempDouble, dataList.get(key))) {
					flag = true;
					break;
				}
				
			}

	    	if(Sum(tempDouble) == false && flag == false) {
	    		dataList.put(count, tempDouble);
		    	count = count + 1;
	    	}
	    		    	
	    }
	    
	    return dataList;

	} 

	
	
	

	
	public static boolean Equality(double[] a, double[] b) {
		
		boolean flag = true;
		for(int i = 0; i < a.length; i++) {
			if(a[i] != b[i]) {
				flag = false;
				break;
			}
		}
		
		return flag;
	
	}
	
	

	public static boolean Sum(double[] data) {
		
		for(int i = 0; i < data.length; i++) {
			if(data[i] == 1.0) {
				return false;
			}
		}
		
		return true;
	}
	
	
	public static void main(String[] args) throws Exception {
		
		int n = 20000;
		int l = 100;
		
		try {
			HashMap<Integer, double[]> dataList = ReadDataJeterWithoutDup(n, l);
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		
		HashMap<Integer, double[]> dataList = ReadDataJeterWithoutDup(n, l);
		HashMap<Integer, ArrayList<Integer>> setHashMap = DataOutsource.vectorMapToSetList(dataList);
		
		ArrayList<Integer> set = new ArrayList<Integer>();
		
		for (Integer key: dataList.keySet()) {
			set.add(setHashMap.get(key).size());
//			System.out.println(setHashMap.get(key).size());
		}
		
		System.out.println("finish");
		Collections.sort(set);
		
		System.out.println(set.get(0));
		System.out.println(set.get(set.size() - 1));
		
		
	}
	
}
