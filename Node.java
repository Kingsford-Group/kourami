public class Node{
        
    public Node(char b, int ci){
	this.base = Character.toUpperCase(b);
	this.colIndex = ci;
	this.iBase = this.base2Index();
    }

    public Node(int ib, int ci){
	
    }
    
    public Node(Base b){
	this(b.getBase(), b.getColPos());
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
	return -1;
    }

    public String toString(){
	return "[" + base + "," + colIndex + "]";
    }
    
    public int getNumPath(){
	return this.numPath;
    }

    public void setNumPath(int n){
	this.numPath = n;
    }
    
    public void updateNumPath(Node pre){
	this.numPath += pre.getNumPath();
    }

    private char base;
    private int iBase;
    private int colIndex;
    
    private int numPath;
}
