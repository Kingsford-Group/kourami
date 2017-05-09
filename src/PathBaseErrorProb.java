/*
Part of Kourami HLA typer/assembler
(c) 2017 by  Heewook Lee, Carl Kingsford, and Carnegie Mellon University.
See LICENSE for licensing.
*/
public class PathBaseErrorProb{

    private char[] pathBases;    
    //public char[][] readBases;
    private double[][] errorProb;
    

    public int numReads(){
	return this.errorProb.length;
    }

    public double[] getNthReadErrorProb(int n){
	return this.errorProb[n];
    }

    public char[] getBases(){
	return this.pathBases;
    }
    
    public double[][] getErrorProb(){
	return this.errorProb;
    }

    public PathBaseErrorProb(int numReads, int len){
	this.pathBases = new char[len];
	//this.readBases = new char[numReads][len];
	this.errorProb = new double[numReads][len];
    }
    
    public void addPathBases(char b, int pos){
	this.pathBases[pos] = b;
    }

    public void add(double e, int readIndex, int pos){
	//this.pathBases[readIndex][pos] = b;
	//this.readBases[readIndex][pos] = b;
	this.errorProb[readIndex][pos] = e;
    }
    
}
