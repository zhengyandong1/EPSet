package SetSim;

import java.util.ArrayList;
import java.util.HashMap;

public class QueryPlain {

	
	public static HashMap<Integer, ArrayList<Integer>> queryPlain(HashMap<Integer, double[]> dataList, ArrayList<Integer> querySet, int tau1, int tau2) {
		
		HashMap<Integer, ArrayList<Integer>> result = new HashMap<Integer, ArrayList<Integer>>();
		
		HashMap<Integer, ArrayList<Integer>> dataListInterger = DataOutsource.vectorMapToSetList(dataList);
		
		for (Integer key : dataListInterger.keySet()) {
			
			if(Trapdoor.jaccardComp(dataListInterger.get(key), querySet) >= tau1*1.0/tau2) {
				result.put(key, dataListInterger.get(key));
			}
			
		}
		
		return result;
		
	}
	
}
