public class Block{

    private int sPos;// start position of MSA
    private StringBuffer sequence;
    
    public Block(int s, String seq){
	this.sPos = s;
	this.sequence = new StringBuffer(seq);
    }
}

class Base{
    
}

class Base{
    
    public char base;
    public int offset;
    public int pos; //actual base position in this sequence
    public boolean exon;
    public int intrexonNum;
    
    public Base(char b, int o, int p, boolean e, int ien){
	this.base = b;
	this.offset = o;
	this.pos = p;
	this.exon = e;
	this.intrexonNum = ien;
    }
    
    

    public ArrayList<Base> loadBlock(String blockS, int boundaryoffset, boolean isExon, int ieNum){
	int baseoffset = 0;
	ArrayList<Base> tmp = new ArrayList<Base>();
	int baseCounts = 0;
	for(int i=0; i<blockS.length(); i++){
	    char curbase = blockS.charAt(i);
	    if(Base.isBase(curbase))
		tmp.add(new Base(curbase, baseoffset, boundaryoffset + (++baseCounts), isExon, ieNum));
	    else if(curbase == '.'){//deletion
		tmp.add(new Base(curbase, 
	    }
	}
    }

    public static boolean isBase(char base){
	if(base == 'A' || base == 'C' || base == 'G' || base = 'T'
	   || base == 'a' || base == 'c' || base == 'g' || base = 't'
	   ){
	    return true;
	}
	return false;
    }
	
}
