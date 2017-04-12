public class Result{
    public int finaScore;
    public int alignedLen;
    public int s1len;
    public int s2len;
    public int identicalLen;
    public double identity;
    public StringBuffer outputBuffer;
    public String hit;

    public Result(int fs, int al, int s1l, int s2l, int il, double id, StringBuffer sb, String hit){
	this.finaScore = fs;
	this.alignedLen = al;
	this.s1len = s1l;
	this.s2len = s2l;
	this.identicalLen = il;
	this.identity = id;
	this.outputBuffer = sb;
	this.identity = this.identicalLen*1.0/(s1l>=s2l ? s1l : s2l); // updated identity using the lenght of longer sequence as its denominator
	this.hit = hit;
    }
    //perfect result
    public Result(int len, String hit){
	this(len, len, len, len, len, 1.0d, new StringBuffer(""), hit);
    }

    public String getHit(){
	return this.hit;
    }

    public String toAlignmentString(){
	return this.outputBuffer.toString();
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
