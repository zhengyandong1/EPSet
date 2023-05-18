package SetSim;

import java.math.BigInteger;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Random;

import SetSim.KDTree.EncNode;

public class Le {
	
	public static SHE paramters = MainFinal.paramters;
	
	public static BigInteger cipherZero1 = MainFinal.cipherZero1;
	public static BigInteger cipherZero2 = MainFinal.cipherZero2; 
	public static BigInteger cipherMinusOne = MainFinal.cipherMinusOne;
	
	
	
	public class TrapdoorLe{
		
		public BigInteger[] trapdoorLe;
		public BigInteger tau1;
		public BigInteger tau2;
		public int sizeQ;
		
		public TrapdoorLe(BigInteger[] trapdoorLe, BigInteger tau1, BigInteger tau2, int sizeQ) {
			
			this.trapdoorLe = trapdoorLe;
			this.tau1 = tau1;
			this.tau2 = tau2;
			this.sizeQ = sizeQ;
			
		}
		
	}
	
	public class Ciphertext{
		
		public BigInteger[] encX;
		public int sizeX;
		
		public Ciphertext(BigInteger[] encX, int sizeX) {
			
			this.encX = encX;
			this.sizeX = sizeX;
		}
		
	}
	
	public static Le le = new Le();
	
	public static ArrayList<Ciphertext> encDataLe(int n, int overElementSize) throws Exception {
		
		int l = 100;
		
		Random random = new Random();
		
		HashMap<Integer, double[]> dataList = Data.ReadDataJeterWithoutDup(n, l);
		
		HashMap<Integer, int[]> dataListIntArray = new HashMap<Integer, int[]>();
		
		for (Integer key : dataList.keySet()) {
			
			double[] data = dataList.get(key);
			
			int[] dataIntArray = new int[overElementSize];
			for(int i = 0; i < dataIntArray.length; i++) {
				if(i < l) {
					
					dataIntArray[i] = (int) data[i];
				}else {
					dataIntArray[i] = 0;
				}
				
			}
			
			dataListIntArray.put(key, dataIntArray);
 			
		}
		
		ArrayList<Ciphertext> encDataSet = new ArrayList<Ciphertext>();
		
		for (Integer key : dataListIntArray.keySet()) {
			
			BigInteger[] encRecord = SHE.encVector(dataListIntArray.get(key), paramters);
			int xSize = random.nextInt();
			Ciphertext cipher = le.new Ciphertext(encRecord, xSize);
			encDataSet.add(cipher);
			
		}
		
		return encDataSet;
		
	}
	
	
	public static TrapdoorLe trapdoorGen(int overElementSize) throws Exception {
		
		int tau1 = 90;
		int tau2 = 100;
		
		Random rnd = new Random();
		
		int[] qSet = new int[overElementSize];
		BigInteger[] encQSet = new BigInteger[overElementSize];
		
		
		for(int i = 0; i < overElementSize; i++) {
			qSet[i] = rnd.nextInt();
			encQSet[i] = SHE.encByCiphertexts(qSet[i], cipherZero1, cipherZero2, paramters);
		}
		
		BigInteger encTau1 = SHE.encByCiphertexts(tau1, cipherZero1, cipherZero2, paramters);
		BigInteger encTau2 = SHE.encByCiphertexts(tau2, cipherZero1, cipherZero2, paramters);
		
		int qSetSize = rnd.nextInt();
		
		TrapdoorLe trapdoorLe = le.new TrapdoorLe(encQSet, encTau1, encTau2, qSetSize);
		
		return trapdoorLe;
		
	}
	
	
	public static double trapdoorGenEvalOne(int overElementSize) throws Exception {
		
		int tau1 = 90;
		int tau2 = 100;
		
		Random rnd = new Random();
		
		int[] qSet = new int[overElementSize];
		BigInteger[] encQSet = new BigInteger[overElementSize];
		
		double t = 0;
		int cycle = 1000;
		
		for(int u = 0; u < cycle; u++) {
			double t1 = System.nanoTime();
			for(int i = 0; i < overElementSize; i++) {
				qSet[i] = rnd.nextInt();
				encQSet[i] = SHE.encByCiphertexts(qSet[i], cipherZero1, cipherZero2, paramters);
			}
			
			BigInteger encTau1 = SHE.encByCiphertexts(tau1, cipherZero1, cipherZero2, paramters);
			BigInteger encTau2 = SHE.encByCiphertexts(tau2, cipherZero1, cipherZero2, paramters);
			
			
			int qSetSize = rnd.nextInt();
			
			TrapdoorLe trapdoorLe = le.new TrapdoorLe(encQSet, encTau1, encTau2, qSetSize);
			double t2 = System.nanoTime();
			
			if(u >= 0.8*cycle) {
				t = t + (t2 - t1);
			}

		}
		
		
		
		return t*5/cycle/1000000;
		
	}
	
	
	public static double queryLeEvalOne(int n, int Uset) throws Exception {
		
		BigInteger N = paramters.sheN;
		Random random = new Random();
		
		double t = 0;
		int cycle = 1000;
		
		for(int i = 0; i < cycle; i++) {
			
			int[] record = new int[Uset];
			for(int j = 0; j < Uset; j++) {
				record[j] = random.nextInt();
			}
			BigInteger[] encRecord = SHE.encVector(record, paramters);
			TrapdoorLe encQset = trapdoorGen(Uset);
			BigInteger encTau1 = SHE.enc(random.nextInt(), paramters);
			BigInteger encTau2 = SHE.enc(random.nextInt(), paramters);
			int SizeX = random.nextInt();
			int SizeQ = random.nextInt();
			
			double t1 = System.currentTimeMillis();
			// Computation at S1
			BigInteger sum = innerProduct(encRecord, encQset.trapdoorLe);
			BigInteger numer = sum.multiply(encTau2);
			BigInteger denumer = BigInteger.valueOf(SizeX + SizeQ);
			denumer = denumer.add(cipherMinusOne.multiply(sum).mod(N)).mod(N);
			denumer = denumer.multiply(encTau2).mod(N);
			
			BigInteger result = numer.add(cipherMinusOne.multiply(denumer).mod(N)).mod(N);
			
			result = SHE.addRandom(result, paramters);
			
			// Computation at S2
			BigInteger val = SHE.dec(result, paramters);
			double t2 = System.currentTimeMillis();
			
			if(i >= cycle*0.8) {
				t = t + (t2 - t1);
			}
		}
		
		t = n*t*5/cycle;
		
		return t;
		
	}
	
	
	public static BigInteger innerProduct(BigInteger[] encRecord, BigInteger[] encSet) {
		
		BigInteger N = paramters.sheN;
		
		BigInteger sum = BigInteger.ZERO;
		for(int i = 0; i < encRecord.length; i++) {
			
			sum = sum.add(encRecord[i].multiply(encSet[i]).mod(N)).mod(N);
			
		}
		
		return sum;
	}
	
	public static double encDatEval(int n, int overElementSize) throws Exception {
		
		double t = 0;
		
		int l = 100;
		
		Random random = new Random();
		
		int cycle = 1000;
		for(int i = 0; i < cycle; i++) {
			
			int[] record = new int[overElementSize];
			for(int j = 0; j < overElementSize; j++) {
				record[j] = random.nextInt();
			}
			double t1 = System.currentTimeMillis();
			BigInteger[] encRecord = SHE.encVector(record, paramters);
			int xSize = random.nextInt();
			Ciphertext cipher = le.new Ciphertext(encRecord, xSize);
			double t2 = System.currentTimeMillis();
			
			if(i >= cycle*0.8) {
				t = t + (t2 - t1);
			}
		}
	
		return n*t*5/cycle;
		
	}
	
	public static void encDataMultiEval() throws Exception {
		
		int n = 20000;
		int[] overElementSize = {200, 400, 600, 800, 1000};
		
		
		for (int i = 0; i < overElementSize.length; i++) {
			
			double t = encDatEval(n, overElementSize[i]);
			System.out.printf("overElementSize = %d, time = %f\n", overElementSize[i], t);
		}
		
		int[] nArray = {4000, 8000, 12000, 16000,  20000};;
		int overElementSizeVal = 200;
		
		for (int i = 0; i < nArray.length; i++) {
			
			double t = encDatEval(nArray[i], overElementSizeVal);
			System.out.printf("n = %d, time = %f\n", nArray[i], t);
		}
		
	}
	
	
	public static void trapdoorGenEval() throws Exception {
		
		
		int[] overElementSize = {200, 400, 600, 800, 1000};
		
		for (int i = 0; i < overElementSize.length; i++) {
			
			double t = trapdoorGenEvalOne(overElementSize[i]); 
			System.out.printf("overElementSize = %d, time = %f\n", overElementSize[i], t);
		}	

	}

	public static void queryLeEval() throws Exception {
		
		
		int[] n = {4000, 8000, 12000, 16000, 20000}; 
		int[] overElementSize = {200, 400, 600, 800, 1000};
		
		for (int i = 0; i < overElementSize.length; i++) {
			System.out.printf("overElementSize = %d\t", overElementSize[i]);
			for(int j = 0; j < n.length; j++) {
				double t = queryLeEvalOne(n[j], overElementSize[i]);
				System.out.printf("%f\t", t);
			}
			System.out.println();
			
		}	

	}
	

}
