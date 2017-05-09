/*
Part of Kourami HLA typer/assembler
(c) 2017 by  Heewook Lee, Carl Kingsford, and Carnegie Mellon University.
See LICENSE for licensing.
*/
import it.unimi.dsi.fastutil.ints.Int2IntMap;
import it.unimi.dsi.fastutil.ints.Int2IntOpenHashMap;
import it.unimi.dsi.fastutil.objects.ObjectIterator;
import it.unimi.dsi.fastutil.ints.IntIterator;

public class CustomHashMap extends Int2IntOpenHashMap{
    
    
    public CustomHashMap(){
	super();
	this.defaultReturnValue(-1);//default phred value of -1 to indicate no entry in the hash.
    }
    
    public CustomHashMap(int expected){
	super(expected);
    }

    public CustomHashMap(int[] k, int[] v){
	super(k, v);
    }
    
    public CustomHashMap(int[] k, int[] v, float f){
	super(k, v, f);
    }

    public CustomHashMap(Int2IntMap m){
	super(m);
    }
    
    public CustomHashMap(Int2IntMap m, float f){
	super(m, f);
    }

    public CustomHashMap(int expected, float f){
	super(expected, f);
    }

    public void printReads(){
	IntIterator itr = this.keySet().iterator();
	HLA.log.append("{");
	while(itr.hasNext())
	    HLA.log.append(" (" +itr.nextInt() + ") ");
	HLA.log.appendln("}");

    }
    
    /* performs union of two HashMap based on Keys */
    public boolean union(CustomHashMap other){
	boolean modified = false;
	for(Int2IntMap.Entry otherEntry : other.int2IntEntrySet()){
	    if(!this.containsKey(otherEntry.getIntKey())){
		this.put(otherEntry.getIntKey(), otherEntry.getIntValue());
		modified = true;
	    }
	}
	return modified;
    }
    
    /* performs intersection of two HashMap based on Keys */
    public boolean intersection(CustomHashMap other){
	boolean modified = false;
	IntIterator itr = this.keySet().iterator();
	while(itr.hasNext()){
	    int curInt = itr.nextInt();
	    if(!other.containsKey(curInt)){
		itr.remove();
		modified = false;
	    }
	}
	return modified;
    }

    /* 
     * Performs a special operation for paired-end
     * Keys of pairing reads are additive inverse to each other.
     * Elements are kepts as long as one additive inverse is in each set.
     *
     */
    public boolean intersectionPE(CustomHashMap other){
	boolean modified = false;
	
	IntIterator itr = this.keySet().iterator();
	while(itr.hasNext()){
	    int key = itr.nextInt();
	    int keyInv = 0 - key;
	    if(!other.containsKey(key) && !other.containsKey(keyInv)){
		itr.remove();
		modified = true;
	    }
	}

	for(int otherKey : other.keySet()){
	    int otherKeyInv = 0-otherKey;
	    boolean intst = this.containsKey(otherKey); //does this contain otherKey?
	    boolean intstPrime = this.containsKey(otherKeyInv); //does this contain -(otherkey)?
	    
	    if(!intst && intstPrime){
		this.put(otherKey, other.get(otherKey));
		modified = true;
	    }
	}
	return modified;
    }

    
    //union of this and (intersectionPE of this and other)
    //this just add PE reads of other that is missing in this readset
    public boolean addPEReads(CustomHashMap other){
	boolean modified = false;
	for(int otherKey : other.keySet()){
	    int otherKeyInv = 0-otherKey;
	    boolean intst = this.containsKey(otherKey);
	    boolean intstPrime = this.containsKey(otherKeyInv);
	    if(!intst && intstPrime){
		this.put(otherKey, other.get(otherKey));
		modified = true;
	    }
	}
	return modified;
    }


    public boolean removeAll(CustomHashMap other){
	boolean modified = false;
	IntIterator itr = other.keySet().iterator();
	int curKey;
	while(itr.hasNext()){
	    curKey = itr.next();
	    if(this.containsKey(curKey)){
		this.remove(curKey);
		modified = true;
	    }
	}
	return modified;
    }

    public boolean addAll(CustomHashMap other){
	return this.union(other);
    }
    
    public boolean retainAll(CustomHashMap other){
	return this.intersection(other);
    }

    /*    public CustomHashMap deepCopy(){
	CustomHashMap tmp = new CustomHashMap();
	for(Int2IntMap.Entry thisEntry : this.int2IntEntrySet())
	    tmp.put(thisEntry.getIntKey(), thisEntry.getIntValue());
	    return tmp;
	}*/
    
    public CustomHashMap clone(){
	return (CustomHashMap)super.clone();
    }

    public void printKeys(){
	IntIterator itr = this.keySet().iterator();
	HLA.log.append("{");
	if(this.size() > 0)
	    HLA.log.append(itr.nextInt());
	while(itr.hasNext())
	    HLA.log.append("," + itr.nextInt());
	HLA.log.appendln("}");
    }
}
