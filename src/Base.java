
public class Base{
    
    public char base; //AaCcGgTT.
    public int iBase; //0 1 2 3 4 A C G T -, 5 for others.
    public int basePos; //1-based basePosition
    public int colPos;  //1-based columnPosition
    public int base2colOffset; 
    public boolean exon;
    public int intronExonNumber; //1-based --> if exon, and intronExonNumber is 1 : Exon1
    public int frame;

    public Base deepCopy(){
	return new Base(this.base, this.iBase, this.basePos, this.colPos, this.base2colOffset, this.exon, this.intronExonNumber, this.frame);
    }

    public Base deepCopyWithOffset(int offsetDifference){
	return new Base(this.base, this.iBase, this.basePos - offsetDifference, this.colPos, this.base2colOffset + offsetDifference, this.exon, this.intronExonNumber, this.frame);
    }

    public String toString(){
	return "b:" + this.base + "|ib:" + iBase + "|bp:" + basePos + "|cp:" + colPos + "|b2co:" + base2colOffset + "|exon?:" + exon + "|ieNum:" + intronExonNumber + "|frame:" + frame;
    }

    public Base(char b, int i, int bp, int cp, int b2co, boolean e, int ien, int f){
	this.base = b;
	this.iBase = i;
	this.basePos = bp;
	this.colPos = cp;
	this.base2colOffset = b2co;
	this.exon = e;
	this.intronExonNumber = ien;
	this.frame = f;
    }

    public Base(char b, int bp, int cp, int b2co, boolean e, int ien){
	this.base = b;
	this.iBase = Base.char2ibase(b);
	this.basePos = bp;
	this.colPos = cp;
	this.base2colOffset = b2co;
	this.exon = e;
	this.intronExonNumber = ien;
	this.frame = -1;

	//if(this.basePos != (this.colPos - this.base2colOffset) ){
	//   System.err.println("Coordinates don't match :" + this.toString());
	//}
    }
    
    public static char ibase2char(int i){
	if(i == 0)
	    return 'A';
	else if(i == 1)
	    return 'C';
	else if(i == 2)
	    return 'G';
	else if(i == 3)
	    return 'T';
	else if(i == 4)
	    return '.';
	else
	    return 'N';
    }
    
    public static int char2ibase(char c){
	if(c == 'A' || c == 'a')
	    return 0;
	else if(c == 'C' || c == 'c')
	    return 1;
	else if(c == 'G' || c == 'g')
	    return 2;
	else if(c == 'T' || c == 't')
	    return 3;
	else if(c == '.' || c == '-')
	    return 4;
	else
	    return 5;
    }

    public boolean isMatch(int ib){
	if(ib == this.iBase)
	    return true;
	return false;
    }

    public boolean isMatch(char b){
	if(this.base == b || this.base == Character.toUpperCase(b))
	    return true;
	return false;
    }

    public static boolean isBase(int ib){
	if(ib >-1 && ib<4 || ib == 5)
	    return true;
	return false;
    }

    public boolean isBase(){
	if(this.iBase > -1 && this.iBase < 4 || this.iBase == 5)
	    return true;
	return false;
    }

    public static boolean isBase(char base){
	if(base == 'A' || base == 'C' || base == 'G' || base == 'T'
	   || base == 'a' || base == 'c' || base == 'g' || base == 't'
	   || base == 'N' || base == 'n'){
	    return true;
	}
	return false;
    }
    
    public static boolean isGap(int ib){
	if(ib == 4)
	    return true;
	return false;
    }

    public static boolean isGap(char base){
	if(base == '.')
	    return true;
	return false;
    }

    public char getBase(){
	return this.base;
    }

    public int getIBase(){
	return this.iBase;
    }

    public Character getBaseUpperObj(){
	return Character.valueOf(Character.toUpperCase(this.base));
    }

    public int getBasePos(){
	return this.basePos;
    }
    
    public int getColPos(){
	return this.colPos;
    }
    
    public int base2colOffset(){
	return this.base2colOffset;
    }

    public int getFrame(){
	return this.frame;
    }
    
    public int getIntronExonNumber(){
	return this.intronExonNumber;
    }
    
    public boolean isExon(){
	return this.exon;
    }
}

