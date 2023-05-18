package SetSim;
import java.math.BigInteger;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;

public class Trapdoor {
	
	public static int k = MainFinal.k;
	public static int scale = MainFinal.scale; 
	public static BigInteger cipherZero1 = MainFinal.cipherZero1;
	public static BigInteger cipherZero2 = MainFinal.cipherZero2; 
	
	public class FilterTrapdoor{
		
		public BigInteger[] leftRange;
		public BigInteger[] rightRange;
		
		public FilterTrapdoor(BigInteger[] leftRange, BigInteger[] rightRange) {
			
			this.leftRange = leftRange;
			this.rightRange = rightRange;
		}
		
	}
	
	public class RefinementTrapdoor{
		
		public ArrayList<BigInteger> encQ;
		public BigInteger encTau1;
		public BigInteger encTau2;
		
		
		public RefinementTrapdoor(ArrayList<BigInteger> encQ, BigInteger encTau1, BigInteger encTau2) {
			
			this.encQ = encQ;
			this.encTau1 = encTau1;
			this.encTau2 = encTau2;
		}
	}
	
	public static Trapdoor trapdoor = new Trapdoor();
	
	
	
	
	public static FilterTrapdoor genFilterTrapdoor(ArrayList<Integer> qSet, HashMap<Integer, ArrayList<Integer>> pivots, int tau1, int tau2, SHE paramters) {
		
		double tau = tau1*1.0/tau2;
		k = MainFinal.k;
		
		int[] leftDisVector = new int[k];
		int[] rightDisVector = new int[k];
		
		for (int i = 0; i < pivots.keySet().size(); i++) {
			
			
			double sim = jaccardComp(qSet, pivots.get(i));
			
			double dis = 1 - sim;
			
			leftDisVector[i] = (int)((dis - (1 - tau))*scale);
			rightDisVector[i] = (int) ((dis + (1 - tau))*scale);
			
		}
		
		
		BigInteger[] leftRange = SHE.encVectorByCiphertexts(leftDisVector, cipherZero1, cipherZero2, paramters);
		BigInteger[] rightRange = SHE.encVectorByCiphertexts(rightDisVector, cipherZero1, cipherZero2, paramters);
		
		FilterTrapdoor filterTrapdoor = trapdoor.new FilterTrapdoor(leftRange, rightRange);
		
		return filterTrapdoor;
		
	}
	
	public static RefinementTrapdoor genRefinementTrapdoor(ArrayList<Integer> qSet, int tau1, int tau2, SHE paramters) {
		
		ArrayList<BigInteger> encQ = SHE.encListByCiphertexts(qSet, cipherZero1, cipherZero2, paramters);
		BigInteger encTau1 = SHE.encByCiphertexts(tau1, cipherZero1, cipherZero2, paramters);
		BigInteger encTau2 = SHE.encByCiphertexts(tau2, cipherZero1, cipherZero2, paramters);
		
		RefinementTrapdoor refinementTrapdoor = trapdoor.new RefinementTrapdoor(encQ, encTau1, encTau2);
		
		return refinementTrapdoor;
		
	}
	
	public static double jaccardComp(ArrayList<Integer> qSet, ArrayList<Integer> xSet) {
		
		int sum = 0;
		
		for(int i = 0; i < qSet.size(); i++) {
			for(int j = 0; j < xSet.size(); j++) {
				if(qSet.get(i) == xSet.get(j)) {
					sum = sum + 1;
					break;
				}
			}
		}
		
		return sum*1.0/(qSet.size() + xSet.size() - sum);
		
	}
	
	public static void main(String[] args) {
		
		
		ArrayList<Integer> qSet = new ArrayList<>();
		
		qSet.add(1);
		qSet.add(2);
		qSet.add(3);
		qSet.add(5);
		qSet.add(7);
		
		ArrayList<Integer> xSet = new ArrayList<>();
		
		xSet.add(4);
		xSet.add(1);
		xSet.add(5);
		xSet.add(3);
		xSet.add(9);
		
		System.out.println(jaccardComp(qSet, xSet));
		
		
	}
	

}
