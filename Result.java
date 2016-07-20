public class Result{
    public int finaScore;
    public int alignedLen;
    public int s1len;
    public int s2len;
    public int identicalLen;
    public double identity;
    
    public Result(int fs, int al, int s1l, int s2l, int il, double id){
	this.finaScore = fs;
	this.alignedLen = al;
	this.s1len = s1l;
	this.s2len = s2l;
	this.identicalLen = il;
	this.identity = id;
    }

    public String toString(){
	return finaScore + "\t" + alignedLen + "\t"+ s1len + "\t"+ s2len + "\t"+ identicalLen + "\t"+ identity;
    }

    public double getIdentity(){
	return this.identity;
    }
    
    public int getIdenticalLen(){
	return this.identicalLen;
    }
    
    public int getAlignedLen(){
	return this.alignedLen;
    }
    
    public int getScore(){
	return this.finaScore;
    }
}
