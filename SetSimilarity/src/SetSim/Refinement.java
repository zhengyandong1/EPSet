package SetSim;

import java.math.BigInteger;
import java.util.ArrayList;
import java.util.Random;
import java.util.concurrent.ForkJoinPool;

import SetSim.Trapdoor.RefinementTrapdoor;

public class Refinement {
	
	public static SHE paramters = MainFinal.paramters;
	
	public static double leafnodeQuery(int sizeSet1, int sizeSet2) {
		
		double t = 0;
		int cycle = 1000;
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
	
	
	
}
