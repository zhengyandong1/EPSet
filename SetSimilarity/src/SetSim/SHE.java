package SetSim;
import java.awt.MultipleGradientPaint.CycleMethod;
import java.math.BigInteger;
import java.security.SecureRandom;
import java.util.ArrayList;
import java.util.Arrays;

import com.sun.security.auth.NTDomainPrincipal;


public class SHE {
	
	
	// paramters of SHE technique
	
	public int k0; // bit length of p and q
	public int k1; // bit length of plaintext space
	public int k2; // bit length of L
	
	public BigInteger shep;		
	public BigInteger sheq;
	public BigInteger sheN;
	public BigInteger sheL;
	public static SecureRandom rnd;
	
	
	public SHE(int k0, int k1, int k2, BigInteger shep, BigInteger sheq, BigInteger sheN, BigInteger sheL) {
		
		this.k0 = k0;
		this.k1 = k1;
		this.k2 = k2;
		this.shep = shep;
		this.sheq = sheq;
		this.sheN = sheN;
		this.sheL = sheL;
	}
	
	// generate the public and secret keys of the system 
	public static SHE keyGen(int k0, int k1, int k2) {
		
		rnd = new SecureRandom();
		
		BigInteger shep = new BigInteger(k0, 40, rnd); // Certainty = 40		
		BigInteger sheq = new BigInteger(k0, 40, rnd); // Certainty = 40
		BigInteger sheN = shep.multiply(sheq);
		BigInteger sheL = new BigInteger(k2, rnd); 
		
		SHE paramters = new SHE(k0, k1, k2, shep, sheq, sheN, sheL);
		
		return paramters;
		
	}
	
	// encrypt an integer value
    public static BigInteger enc(int val, SHE paramters) {
    	
    	int k0 = paramters.k0; // bit length of p and q
    	int k2 = paramters.k2; // bit length of L
    	
    	BigInteger p = paramters.shep;
    	BigInteger L = paramters.sheL;
    	BigInteger N = paramters.sheN;
    	
    	BigInteger message = BigInteger.valueOf(val);
    	
    	BigInteger r2 = new BigInteger(k2, rnd);  // Certainty = 40
    	BigInteger r0 = new BigInteger(k0, rnd);  // Certainty = 40
    	
    	BigInteger c = message.add(r2.multiply(L));
    	c = c.multiply(BigInteger.ONE.add(r0.multiply(p))).mod(N);
    	
    	return c;
    	
    }
    
    // decrypt an integer value
    public static BigInteger dec(BigInteger c, SHE paramters) {
    	
    	BigInteger p = paramters.shep;
    	BigInteger L = paramters.sheL;
    	
    	BigInteger message = c.mod(p).mod(L);
    	
    	if(message.compareTo(L.divideAndRemainder(BigInteger.TWO)[0]) == 1) {
    		
    		message = message.subtract(L);
 
    	}
    	
    	return message;
    	
    }     
    
    // encrypt an integer vector
    public static BigInteger[] encVector(int[] vector, SHE paramters) {
    	
    	int d = vector.length;
    	
    	BigInteger[] cVector = new BigInteger[d];
    	
    	for(int i = 0; i < d; i++) {
    		cVector[i] = enc(vector[i], paramters);
    	}
    	
    	return cVector;
		
	}
    
    // decrypt an integer vector
    public static BigInteger[] decVector(BigInteger[] cVector, SHE paramters) {
     	
     	int d = cVector.length;
     	
     	BigInteger[] vector = new BigInteger[d];
     	
     	for(int i = 0; i < d; i++) {
     		vector[i] = dec(cVector[i], paramters);
     	}
     	
     	return vector;
 		
 	}
    
    
    // encrypt a list of integers
    public static ArrayList<BigInteger> encIntegerList(ArrayList<Integer> list, SHE paramters) {
    	
    	int d = list.size(); 
    	
    	ArrayList<BigInteger> encRecordList = new ArrayList<BigInteger>();
    	
    	for(int i = 0; i < d; i++) {
    		
    		encRecordList.add(enc(list.get(i), paramters));
    		
    	}
    	
    	return encRecordList;
		
	}
    
    public static ArrayList<BigInteger> decIntegerList(ArrayList<BigInteger> cipherList, SHE paramters) {
    	
    	int d = cipherList.size(); 
    	
    	ArrayList<BigInteger> decList = new ArrayList<BigInteger>();
    	
    	for(int i = 0; i < d; i++) {
    		
    		decList.add(dec(cipherList.get(i), paramters));
    		
    	}
    	
    	return decList;
	}

    
    
    // encrypt an integer by ciphertexts
    public static BigInteger encByCiphertexts(int val, BigInteger cipherZero1, BigInteger cipherZero2, SHE paramters) {
    	
    	int k2 = paramters.k2;
    	BigInteger N = paramters.sheN;
    	
    	BigInteger r1 = new BigInteger(k2, rnd);  // Certainty = 40
    	BigInteger r2 = new BigInteger(k2, rnd);  // Certainty = 40 
    	
    	BigInteger c = cipherZero1.multiply(r1).mod(N);
    	c = c.add(cipherZero2.multiply(r2).mod(N));
    	c = c.add(BigInteger.valueOf(val)).mod(N);
    	
    	return c;
    	
    }
    
    // encrypt a vector by ciphertexts
    public static BigInteger[] encVectorByCiphertexts(int[] vector, BigInteger cipherZero1, BigInteger cipherZero2, SHE paramters) {
    	
    	int d = vector.length;
    	BigInteger[] cVector = new BigInteger[d];
    	
    	for(int i = 0; i < d; i++) {
    		cVector[i] = encByCiphertexts(vector[i], cipherZero1, cipherZero2, paramters);
    	}
    	
    	
    	return cVector;
    	
    }    
    
    // encrypt a list by ciphertexts
    public static ArrayList<BigInteger> encListByCiphertexts(ArrayList<Integer> dataList, BigInteger cipherZero1, BigInteger cipherZero2, SHE paramters) {
    	
    	int d = dataList.size();
    	
    	ArrayList<BigInteger> encDataList = new ArrayList<BigInteger>();
    	
    	for(int i = 0; i < d; i++) {
    		
    		BigInteger ciphertext = encByCiphertexts(dataList.get(i), cipherZero1, cipherZero2, paramters);
    		encDataList.add(ciphertext);
    		
    	}
    	
    	return encDataList;
    	
    }
    
    
    public static BigInteger addRandom(BigInteger c, SHE paramters) {
    	
    	int k1 = paramters.k1;
    	BigInteger N = paramters.sheN;
    	
    	BigInteger r1;
    	BigInteger r2;
    	
    	while(true) {
    		
    		r1 = new BigInteger(k1, rnd);  // Certainty = 40
        	r2 = new BigInteger(k1, rnd);  // Certainty = 40 
        	
        	if(r1.compareTo(r2) == 1) {
        		break;
        	}
    	}
    	
    	c = r1.multiply(c).add(r2).mod(N); 
    	
    	return c;
    	
    }
    
    public static ArrayList<BigInteger> addList(ArrayList<BigInteger> list1, ArrayList<BigInteger> list2, SHE paramters) {
    	
    	BigInteger N = paramters.sheN;
    	
    	ArrayList<BigInteger> list = new ArrayList<BigInteger>();
    	
    	for(int i = 0; i < list1.size(); i++) {
    		list.add(list1.get(i).add(list2.get(i)).mod(N));
    	}
    	
    	return list;
    }
    
    public static ArrayList<Integer> arrayToList(int[] m) {
    	
    	ArrayList<Integer> list = new ArrayList<Integer>();
    	
    	for(int i = 0; i < m.length; i++) {
    		list.add(m[i]);
    	}
    	
    	return list;
    	
    } 

    
    public static void main(String[] args) {
		
    	int k0 = 1024;
    	int k1 = 40;
    	int k2 = 160;
    	
    	// key generation of our scheme
    	SHE paramters = keyGen(k0, k1, k2);
    	
    	
    	// test encrypt val 
    	
//    	int m = 1;
//    	BigInteger c = enc(m, paramters);
//    	System.out.println(dec(c, paramters));
//    	
//    	m = -1;
//    	c = enc(m, paramters);
//    	System.out.println(dec(c, paramters));
    	
    	
//    	// encryption vector
//    	int[] m = {1,-2,3,-4,5,10, -1000, 10000};
//    	BigInteger[] c = encVector(m, paramters);
//    	
//    	for(int i = 0; i < m.length; i++) {
//    		System.out.printf("%d,", dec(c[i], paramters));
//    	}
//    	System.out.println();
    	
//    	// encryption list
//    	int[] m = {1,-2,3,-4,5,10, -1000, 10000};
//    	ArrayList<Integer> mList = arrayToList(m);
//    	
//    	ArrayList<BigInteger> plainlist = encIntegerList(mList, paramters);
//    			
//    	
//    	for(int i = 0; i < m.length; i++) {
//    		System.out.printf("%d,", dec(plainlist.get(i), paramters));
//    	}
//    	System.out.println();
    	
    	
    	
    	
    	
    	
     	
	}

	
}
