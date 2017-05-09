/*
Part of Kourami HLA typer/assembler
(c) 2017 by  Heewook Lee, Carl Kingsford, and Carnegie Mellon University.
See LICENSE for licensing.
*/
public class Result{
    public int finaScore;
    public int alignedLen;
    public int s1len;
    public int s2len;
    public int identicalLen;
    public double identity;
    public StringBuffer outputBuffer;
    public String hit;
    public String gGroupName;
    public boolean perfectMatch;

    public Result(int fs, int al, int s1l, int s2l, int il, double id, StringBuffer sb, String hit, String ggn, boolean ip){
	this.finaScore = fs;
	this.alignedLen = al;
	this.s1len = s1l;
	this.s2len = s2l;
	this.identicalLen = il;
	this.identity = id;
	this.outputBuffer = sb;
	this.identity = this.identicalLen*1.0/(s1l>=s2l ? s1l : s2l); // updated identity using the lenght of longer sequence as its denominator
	this.hit = hit;
	this.gGroupName = ggn;
	this.perfectMatch = ip;
    }
    //imperfect match
    public Result(int fs, int al, int s1l, int s2l, int il, double id, StringBuffer sb, String hit, String ggn){
	this(fs, al, s1l, s2l, il, id
	     , sb, hit, ggn, false);
    }
    //perfect result
    public Result(int len, String hit, String ggn){
	this(len, len, len, len, len, 1.0d
	     , new StringBuffer(""), hit, ggn, true);
    }
    //perfect result
    public Result(int len, HLASequence hs2){
	this(len, len, len, len, len, 1.0d
	     , new StringBuffer(""), hs2.getSequence(), hs2.getGroup().getGroupString(), true);
    }

    public double getPairIdentity(Result other){
	return (this.identicalLen + other.getIdenticalLen())*1.0d / (this.getMaxLen() + other.getMaxLen())*1.0d;
    }

    public int getMaxLen(){
	if(this.s1len >= this.s2len)
	    return this.s1len;
	return this.s2len;
    }

    

    public String getGGroupName(){
	return this.gGroupName;
    }

    public boolean isPerfect(){
	return this.isPerfect();
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
