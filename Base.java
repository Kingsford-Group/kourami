
public class Base{
    
    public char base; //AaCcGgTT.
    public int basePos; //1-based basePosition
    public int colPos;  //1-based columnPosition
    public int base2colOffset; 
    public boolean exon;
    public int intronExonNumber; //1-based --> if exon, and intronExonNumber is 1 : Exon1
    
    public Base(char b, int bp, int cp, int b2co, boolean e, int ien){
	this.base = b;
	this.basePos = bp;
	this.colPos = cp;
	this.base2colOffset = b2co;
	this.exon = e;
	this.intronExonNumber = ien;
    }
    
    public boolean isMatch(char b){
	if(this.base == b || this.base == Character.toUpperCase(b))
	    return true;
	return false;
    }

    public static boolean isBase(char base){
	if(base == 'A' || base == 'C' || base == 'G' || base = 'T'
	   || base == 'a' || base == 'c' || base == 'g' || base = 't'
	   ){
	    return true;
	}
	return false;
    }
    
    public static boolean isGap(char base){
	if(base == '.')
	    return true;
	return false;
    }
    
}

