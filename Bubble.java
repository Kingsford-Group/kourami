import java.io.*;
import java.util.*;

public class Bubble{

    private HLAGraph g;
    //private int start;
    //private int end;
    private ArrayList<Integer> start;
    private ArrayList<Integer> end;
    
    private ArrayList<Node> sNodes;
    private ArrayList<Node> tNodes;
    //    private Node s;
    //private Node t;
    
    private ArrayList<Path> paths;

    public void printBubbleSequenceSizes(){
	for(Path p : this.paths)
	    System.err.print(p.getBubbleSequences().size() + "\t");
	System.err.println();
    }

    public ArrayList<Integer> getStart(){
	return this.start;
    }
    
    public ArrayList<Integer> getEnd(){
	return this.end;
    }
    
    public ArrayList<Node> getSNodes(){
	return this.sNodes;
    }

    public ArrayList<Node> getTNodes(){
	return this.tNodes;
    }


    public ArrayList<Path> getPaths(){
	return this.paths;
    }

    public void printBubbleSequence(){
	for(Path p : this.paths){
	    p.printBubbleSequence(g.getGraph());
	}
    }

    public void initBubbleSequences(){
	for(Path p : this.paths){
	    p.initBubbleSequence(g.getGraph());
	}
    }

    public void printResults(ArrayList<StringBuffer> interBubbleSequences){
	System.err.println("Printing\t" + this.paths.size() + "\tpossible sequences"  );
	for(int i=0; i<this.paths.size(); i++){
	    Path p = this.paths.get(i);
	    ArrayList<StringBuffer> bubbleSequences = p.getBubbleSequences();
	    System.out.println(">>>>>>>ITBS\t" + interBubbleSequences.size() + "\tBS\t" + bubbleSequences.size());
	    System.out.print(interBubbleSequences.get(0).toString());
	    for(int j=1; j<interBubbleSequences.size(); j++){
		System.out.print("<<" + bubbleSequences.get(j-1).toString() + ">>");
		System.out.print(interBubbleSequences.get(j).toString());
	    }
	}
    }
    
    

    //print fractured bubbles
    //returns next startIndex
    public int printResults(ArrayList<StringBuffer> interBubbleSequences, int startIndex){
	int tmpStartIndex = startIndex;
	for(int i=0; i<this.paths.size(); i++){
	    Path p = this.paths.get(i);
	    System.out.println("IntersectionScore:\t" + p.getAvgWeightedIntersectionSum());
	    ArrayList<StringBuffer> bubbleSequences = p.getBubbleSequences();
	    //each bubbleSequence is padded by interBubbleSequences
	    //so we print the first interBubbleSequence.
	    tmpStartIndex = startIndex;
	    System.out.print(interBubbleSequences.get(tmpStartIndex));
	    tmpStartIndex++;
	    
	    for(int j=0; j<bubbleSequences.size(); j++){
		System.out.print(" <" + bubbleSequences.get(j) + "> ");//prints the bubble
		System.out.print(interBubbleSequences.get(tmpStartIndex).toString()); //prints the interBubble
		tmpStartIndex++;
	    }
	    System.out.println();
	}
	return tmpStartIndex;
    }

    /*
    public int printResults(ArrayList<StringBuffer> interBubbleSequences, int startIndex){
	int nextStartIndex = 0;
	int tmpStartIndex = startIndex;
	for(int i=0; i<this.paths.size(); i++){
	    Path p = this.paths.get(i);
	    ArrayList<StringBuffer> bubbleSequences = p.getBubbleSequences();
	    
	    if(startIndex == 0){
		System.out.print(interBubbleSequences.get(startIndex));
		tmpStartIndex=1;
	    }
	    nextStartIndex = tmpStartIndex + bubbleSequences.size();
	    for(int j=0; j<bubbleSequences.size(); j++){
		System.out.print(" <" + bubbleSequences.get(j) + "> ");
		System.out.print(interBubbleSequences.get(tmpStartIndex+j).toString());
	    }
	    System.out.println();
	}
	return nextStartIndex;
    }
    */

    public Bubble(HLAGraph hg, Node s, Node t){
	this.g = hg;
	this.sNodes = new ArrayList<Node>();
	this.tNodes = new ArrayList<Node>();
	this.sNodes.add(s);
	this.tNodes.add(t);
	this.start = new ArrayList<Integer>();
	this.end = new ArrayList<Integer>();
	this.start.add(new Integer(s.getColIndex()));
	this.end.add(new Integer(t.getColIndex()));
	this.paths = new ArrayList<Path>();
	this.decompose(s, t);
	this.removeUnsupported();
    }
    
    //find all ST path in the bubble
    //and remove unsupported paths
    public ArrayList<Path> decompose(Node s, Node t){
	System.err.print("Bubble decomposing...\t");
	this.paths =  this.g.findAllSTPath(s, t);
	System.err.print("Found (" + this.paths.size() + ") possible paths.\n");// + "Removed (");
	this.initPathCounters();
	for(Path p : this.paths){
	    //p.printNumActivePaths();
	    p.includePath();
	}
	
	for(Path p : this.paths){
	    p.computeReadSet(g);
	}
	
	/*int numRemoved = this.removeUnsupported();
	System.err.println(numRemoved + ") paths.");
	System.err.println("Total Left:\t"+ this.paths.size());
	if(numRemoved > 0){
	    System.err.println("After removal:");
	    for(Path p : this.paths){
		p.printNumActivePaths();
	    }
	    }*/
	return this.paths;
    }
    
    public void initPathCounters(){
	for(Path p : this.paths){
	    p.initPathCounter();
	}
    }
    
    public void printPaths(){
	int count = 0;
	for(Path p : this.paths){
	    System.err.print("\tP" + count + "\t");
	    p.printPath();
	    count++;
	}
    }
    
    public int removeUnsupported(){
	System.err.println("[Bubble] unsupported path removal...");
	HashSet<Integer> readHash;
	ArrayList<Integer> removalList = new ArrayList<Integer>();
	for(int i=0; i<this.paths.size(); i++){
	    Path p = this.paths.get(i);
	    if(!p.isSupportedPath()){
		removalList.add(new Integer(i));
		System.err.print("Removing\tPath" + i + "\t");
		p.printPath();
	    }
	}
	
	for(int i=removalList.size() - 1; i >= 0; i--){
	    this.paths.get(removalList.get(i).intValue()).excludePath();
	    this.paths.remove(removalList.get(i).intValue());
	}
	System.err.println("Removed (" + removalList.size() + ") paths and\t(" + this.paths.size() + ") left.");
	return removalList.size();
    }

    public int removeUnsupported2(){
	System.err.println("[Bubble] unsupported path removal...");
	HashSet<Integer> readHash;
	ArrayList<Integer> removalList = new ArrayList<Integer>();
	for(int i=0; i<this.paths.size(); i++){//Path p: this.paths){
	    Path p = this.paths.get(i);
	    ArrayList<CustomWeightedEdge> eList = p.getOrderedEdgeList();
	    boolean inited = false;
	    readHash = null;
	    //int numUnique = 0;
	    StringBuffer bf = new StringBuffer();
	    for(CustomWeightedEdge e : eList){
		bf.append("(" + e.getEdgeId() + ")");
		//System.err.print("|" + e.getNumActivePath() + "|");
		//if(e.isUniqueEdge()){ //we only care about uniq edges
		    
		//numUnique++;
		if(!inited){
		    readHash = e.getReadHashSetDeepCopy();
		    inited = true;
		}else{
		    //update with union after checking intersection
		    //if null, we remove this path
		    

		    if(e.unionAfterCheckingIntersection(readHash) == null){
			removalList.add(new Integer(i));
			break;
		    }
		}
		
	    }
	    System.err.println("Removing\tPath" + i+"\t" + bf.toString());
	    //System.err.print("[" + numUnique +  "]");
	}
	
	for(int i=removalList.size() - 1; i >= 0; i--){
	    this.paths.get(removalList.get(i).intValue()).excludePath();
	    this.paths.remove(removalList.get(i).intValue());
	    
	}
	System.err.println("Removed (" + removalList.size() + ") paths and\t(" + this.paths.size() + ") left.");
	return removalList.size();
    }

    /* same as removeUnsupported except this only cares about unique edges */
    public int removeUnsupportedUniqueEdgeOnly(){
	HashSet<Integer> readHash;
	ArrayList<Integer> removalList = new ArrayList<Integer>();
	for(int i=0; i<this.paths.size(); i++){//Path p: this.paths){
	    Path p = this.paths.get(i);
	    ArrayList<CustomWeightedEdge> eList = p.getOrderedEdgeList();
	    boolean inited = false;
	    readHash = null;
	    int numUnique = 0;
	    for(CustomWeightedEdge e : eList){
		System.err.print("|" + e.getNumActivePath() + "|");
		if(e.isUniqueEdge()){ //we only care about uniq edges
		    
		    numUnique++;
		    if(!inited){
			readHash = e.getReadHashSetDeepCopy();
			inited = true;
		    }else{
			//update with union after checking intersection
			//if null, we remove this path
			if(e.unionAfterCheckingIntersection(readHash) == null){
			    removalList.add(new Integer(i));
			    break;
			}
		    }
		}
	    }
	    System.err.print("[" + numUnique +  "]");
	}
	for(int i=removalList.size() - 1; i >= 0; i--){
	    this.paths.get(removalList.get(i).intValue()).excludePath();
	    this.paths.remove(removalList.get(i).intValue());
	    
	}
	return removalList.size();
    }



    /*
    public void removeUnsupported(){
	HashSet<Integer> readHash;
	ArrayList<Integer> removalList = new ArrayList<Integer>();
	for(int i=0; i<this.paths.size(); i++){//Path p: this.paths){
	    Path p = this.paths.get(i);
	    ArrayList<CustomWeightedEdge> eList = p.getOrderedEdgeList();
	    boolean inited = false;
	    for(CustomWeightedEdge e : eList){
		if(e.isUnique()){//we only care about uniq edges
		    if(!inited){
			readHash = e.getReadHashSet();
			inited = true;
		    }else{
			HashSet<Integer> curHash = e.getReadHashSet();
			if(curHash.retainAll(readHash)){
			    if(curHash.isEmpty()
			}

			if(curHash.retainAll(readHash)){
			    if(curHash.isEmpty())
				removalList.add(new Integer(i));
			}
			}
		}
	    }
	}
	    for(int i=removalList.size() - 1; i >= 0; i--){
	    this.paths.remove(removalList.get(i));
	}
    }
    */
    public boolean mergeBubble(Bubble other){
	ArrayList<Path> paths_new = new ArrayList<Path>();
	//boolean[] tpUsed = new boolean[this.paths.size()];
	//boolean[] opUsed = new boolean[other.getPaths().size()];
	
	/* path used counters */
	int[] tpUsed = new int[this.paths.size()];
	int[] opUsed = new int[other.getPaths().size()];

	
	/* print this paths (DEBUGGIN) */
	for(int i=0;i<this.paths.size();i++){
	    Path tp = this.paths.get(i);
	    System.err.print("TP(" + i + ")\t<readNum:" + tp.getReadSetSize() + ">\t");
	    tp.printInfo();
	}
	
	/* print other paths (DEBUGGIN) */
	for(int i=0;i<other.getPaths().size();i++){
	    Path op = other.getPaths().get(i);
	    System.err.print("OP(" + i + ")\t<readNum:" + op.getReadSetSize() + ">\t");
	    op.printInfo();
	}
	
	ArrayList<int[]> phasedList = new ArrayList<int[]>();
	ArrayList<Integer> intersectionSizes = new ArrayList<Integer>();
	int intersectionSizesSum = 0;
	/* check possible paths TP X OP */
	for(int i=0;i<this.paths.size();i++){
	    Path tp = this.paths.get(i);
	    for(int j=0; j<other.getPaths().size(); j++){
		Path op = other.getPaths().get(j);
		System.err.print("TP(" + i + ")\tX\tOP("+j+"):\t");
		int intersectionSize = tp.isPhasedWith(op);
		if(intersectionSize >= Path.MIN_SUPPORT_PHASING){
		    //if(tp.isPhasedWith(op)){
		    //paths_new.add(tp.mergePaths(op));
		    intersectionSizesSum += intersectionSize;
		    intersectionSizes.add(new Integer(intersectionSize));
		    int[] tmp = new int[2];
		    tmp[0] = i;
		    tmp[1] = j;
		    phasedList.add(tmp);
		    tpUsed[i]++;
		    opUsed[j]++;
		}
	    }
	}
	
	if(phasedList.size() == 0){
	    return false;//System.err.println("CANT PHASE!!!!!!");
	}else{
	    System.err.println("TOTAL of " + phasedList.size() + "\tphased paths.");
	
	    
	    //for(int[] ijs : phasedList){
	    for(int i=0;i<phasedList.size();i++){
		int[] ijs = phasedList.get(i);
		int intersectionSize = intersectionSizes.get(i);
		Path tp = this.paths.get(ijs[0]);
		Path op = other.getPaths().get(ijs[1]);
		//if tp and op are used once, merged path between tp and op is the only PATH
		if(tpUsed[ijs[0]] == 1 && opUsed[ijs[1]] == 1){
		    paths_new.add(tp.mergePathsUnique(op));
		}
		//tp is used once  op is used multiple times, we pass tp readset.
		else if(tpUsed[ijs[0]] == 1 && opUsed[ijs[1]] > 1){
		    paths_new.add(tp.mergePath1toMany(op));
		}
		//tp is used multiple time, op is used once, we pass op readset.
		else if(tpUsed[ijs[0]] > 1 && opUsed[ijs[1]] == 1){
		    paths_new.add(tp.mergePathManyto1(op));
		}
		//tp is used multiple time, op is used multiple times, we pass intersection.
		else if(tpUsed[ijs[0]] > 1 && opUsed[ijs[1]] > 1){
		    paths_new.add(tp.mergePathManytoMany(op));
		}else{
		    System.err.println("SOMETHING IS WRONG. [Bubble.java mergeBubble()]");
		    System.exit(-1);
		}
		paths_new.get(paths_new.size() - 1).updateIntersectionSum(intersectionSize, intersectionSizesSum);

	    }
	    System.err.println(paths_new.size() + "\tphased paths in paths_new");
	    
	    /* edge usage update for TP */
	    for(int i=0; i<this.paths.size();i++){
		if(tpUsed[i] == 0)
		    this.paths.get(i).excludePath();
		else 
		    this.paths.get(i).includePathNTimes(tpUsed[i]-1);
	    }
	    
	    /* edge usage update for OP */
	    for(int i=0; i<other.getPaths().size();i++){
		if(opUsed[i] == 0)
		    other.getPaths().get(i).excludePath();
		else
		    other.getPaths().get(i).includePathNTimes(opUsed[i]-1);
	    }
	    
	    if(paths_new.size() > 0){
		this.paths = paths_new;
		this.start.addAll(other.getStart());
		this.end.addAll(other.getEnd());
		this.sNodes.addAll(other.getSNodes());
		this.tNodes.addAll(other.getTNodes());
	    }
	
	    return true;
	}
    }

    public int getNumPaths(){
	return this.paths.size();
    }
    
}
