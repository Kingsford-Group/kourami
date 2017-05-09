/*
Part of Kourami HLA typer/assembler
(c) 2017 by  Heewook Lee, Carl Kingsford, and Carnegie Mellon University.
See LICENSE for licensing.
*/
//import java.util.HashSet;

public class Node{
        
    public Node(char b, int ci){
	this.base = Character.toUpperCase(b);
	this.colIndex = ci;
	this.iBase = this.base2Index();
	//this.rHash = new HashSet<Integer>();
	//moved rHash to CustomWeightedEdge
    }

    public Node(int ib, int ci){
	
    }

    //moved rHash to CustomWeightedEdge
    /*
    public void addAllReadsFrom(HashSet<Integer> otherRHash){
	this.rHash.addAll(otherRHash);
    }

    public HashSet<Integer> getReadHashSet(){
	return this.rHash;
	}*/
    
    public Node(Base b){
	this(b.getBase(), b.getColPos());
    }

    public int getIBase(){
	return this.iBase;
    }
    
    public char getBase(){
	return this.base;
    }

    public Character getBaseObj(){
	return new Character(this.base);
    }
    
    public int getColIndex(){
	return this.colIndex;
    }
    
    public void setColIndex(int ni){
	this.colIndex = ni;
    }

    public boolean equals(Node other){
	if(this.base == other.getBase()){
	    if(this.colIndex == other.getColIndex())
		return true;
	}
	return false;
    }

    public int base2Index(){
	if(this.base == 'A' || this.base == 'a')
	    return 0;
	else if(this.base == 'C' || this.base == 'c')
	    return 1;
	else if(this.base == 'G' || this.base == 'g')
	    return 2;
	else if(this.base == 'T' || this.base == 't')
	    return 3;
	else if(this.base == '-' || this.base == '.')
	    return 4;
	else
	    return 5;
    }

    public String toString(){
	return "[" + base + "," + colIndex + "]";
    }
    
    /*
    public void incrementNumPathInBubbleFwd(int inc){
	this.numPathInBubbleFwd += inc;
    }
    
    public void incrementNumPathInBubbleRev(int inc){
	this.numPathInBubbleFwd += inc;
    }

    public void setNumInBubbleFwd(int n){
	this.numPathInBubbleFwd = n;
    }
    
    public void setNumInBubbleRev(int n){
	this.numPathInBubbleRev = n;
    }
    
    public int getNumInBubbleFwd(){
	return this.numPathInBubbleFwd;
    }
    
    public int getNumInBubbleRev(){
	return this.numPathInBubbleRev;
    }

    public void initBubblePathsCounters(){
	this.numPathInBubbleFwd = 0;
	this.numPathInBubbleRev = 0;
    }
    */
    private char base;
    private int iBase;
    private int colIndex;
    
    private int numPathInBubbleFwd;
    private int numPathInBubbleRev;

    //moved rHash to CustomWeightedEdge
    /*
    public void addRead(int readNum){
	this.rHash.add(new Integer(readNum));
    }
    
    private HashSet<Integer> rHash;
    */
}
