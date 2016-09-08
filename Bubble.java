import java.io.*;
import java.util.*;

public class Bubble{

    private HLAGraph g;
    //private int start;
    //private int end;
    private ArrayList<Integer> start;
    private ArrayList<Integer> end;
    
    //s- and t- nodes of each bubbles that are merged
    private ArrayList<Node> sNodes;
    private ArrayList<Node> tNodes;
    //    private Node s;
    //private Node t;

    private ArrayList<Integer> bubbleLengths; // when bubbles get merged bubbleLengths keep track of lengths of each bubble being merged. The length of this list is always equal to the number of merged bubbles.

    //keeps track of paths through the whole bubble(merged or not-merged)
    private ArrayList<Path> paths;

    private boolean firstBubble;

    public boolean isFirstBubble(){
	return this.firstBubble;
    }

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
    
    public static String stripPadding(String columnSeq){
	StringBuffer build = new StringBuffer();
	for(int i=0; i<columnSeq.length(); i++){
	    char curchar = columnSeq.charAt(i);
	    if(curchar == '.' || curchar == ' ')
		;
	    else
		build.append(curchar);
	}
	return build.toString();
    }

    //print fractured bubbles
    //returns next startIndex
    public int printResults(ArrayList<StringBuffer> interBubbleSequences, int startIndex, ArrayList<DNAString> sequences, String hlagenename, int superbubbleNumber){
	int tmpStartIndex = startIndex;
	for(int i=0; i<this.paths.size(); i++){
	    Path p = this.paths.get(i);
	    int pathnum = i;
	    DNAString curDNA = new DNAString(hlagenename, superbubbleNumber, pathnum
					     , p.getAvgWeightedIntersectionSum()
					     , p.getProbability());
	    //output.append(">" + hlagenename + "_" + superbubbleNumber + "_" + pathnum + "-" + p.getAvgWeightedIntersectionSum() + ":" + p.getProbability() + "\n");
	    System.out.println("IntersectionScore:\t" + p.getAvgWeightedIntersectionSum() + "\t" + p.getProbability());
	    ArrayList<StringBuffer> bubbleSequences = p.getBubbleSequences();
	    //each bubbleSequence is padded by interBubbleSequences
	    //so we print the first interBubbleSequence.
	    tmpStartIndex = startIndex;
	    if(superbubbleNumber == 0 || this.firstBubble){
		System.out.print(interBubbleSequences.get(tmpStartIndex).toString());
		curDNA.append(Bubble.stripPadding(interBubbleSequences.get(tmpStartIndex).toString()));
		//output.append(Bubble.stripPadding(interBubbleSequences.get(tmpStartIndex).toString()));
		tmpStartIndex++;
	    }
	    
	    
	    for(int j=0; j<bubbleSequences.size(); j++){
		System.out.print(" <" + bubbleSequences.get(j) + "> ");//prints the bubble
		curDNA.append(Bubble.stripPadding(bubbleSequences.get(j).toString()));
		//output.append(Bubble.stripPadding(bubbleSequences.get(j).toString()));

		//if path ends with bubble, we dont need to print the last interBubbleSequence
		if(tmpStartIndex < interBubbleSequences.size()){
		    System.out.print(interBubbleSequences.get(tmpStartIndex).toString()); //prints the interBubble
		    curDNA.append(Bubble.stripPadding(interBubbleSequences.get(tmpStartIndex).toString()));
		    //output.append(Bubble.stripPadding(interBubbleSequences.get(tmpStartIndex).toString()));
		}
		
		tmpStartIndex++;
	    }
	    System.out.println();
	    //curDNA.append("\n");
	    //output.append("\n");
	    sequences.add(curDNA);
	}
	return tmpStartIndex;
    }


    public int mergePathsInSuperBubbles(ArrayList<Path> interBubblePaths, int startIndex, ArrayList<Path> resultPaths, String hlagenename, int superbubbleNumber){
	int tmpStartIndex = startIndex;
	for(int i=0; i<this.paths.size(); i++){
	    Path p = this.paths.get(i);
	    int pathnum = i;
	    Path curPath = new Path();
	    curPath.setProbability(p.getProbability());
	    curPath.setWeightedIntersectionSum(p.getWeightedIntersectionSum());
	    curPath.setMergedNums(p.getMergedNums());
	    tmpStartIndex = startIndex;
	    if(superbubbleNumber == 0 || this.firstBubble){
		curPath.appendAllEdges(interBubblePaths.get(tmpStartIndex));
		tmpStartIndex++;
	    }
	    
	    //for each bubble merged in this path
	    int k=0;
	    for(int j=0; j<this.bubbleLengths.size(); j++){
		int bubbleLength = this.bubbleLengths.get(j);
		int limit = k+bubbleLength;
		for(;k<limit;k++){
		    curPath.appendEdge(p.getNthEdge(k));
		}
		if(tmpStartIndex < interBubblePaths.size()){
		    curPath.appendAllEdges(interBubblePaths.get(tmpStartIndex));
		}
		tmpStartIndex++;
	    }
	    resultPaths.add(curPath);
	}
	return tmpStartIndex;
    }

    /*
    public int mergePathsInSuperBubbles(ArrayList<Path> interBubblePaths, int startIndex, ArrayList<Path> paths, String hlagenename, int superbubbleNumber){
    
	int tmpStartIndex = startIndex;
	for(int i=0; i<this.paths.size(); i++){
	    Path p = this.paths.get(i);
	    int pathnum = i;
	    Path curP = new Path();
	    curP.
	    tmpStartIndex = startIndex;
	    if(superbubbleNumber == 0 || this.firstBubble){
		
	    }
	}
    }
    */
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
	this.firstBubble = false;
	this.g = hg;
	this.sNodes = new ArrayList<Node>();
	this.tNodes = new ArrayList<Node>();
	this.sNodes.add(s);
	this.tNodes.add(t);
	this.start = new ArrayList<Integer>();
	this.end = new ArrayList<Integer>();
	if(s == null){
	    System.err.println("[BUBBLE] start node null");
	}
	this.start.add(new Integer(s.getColIndex()));
	this.end.add(new Integer(t.getColIndex()));
	this.paths = new ArrayList<Path>();
	this.decompose(s, t);
	this.removeUnsupported();
    }

    public Bubble(HLAGraph hg, Node s, Node t, boolean firstBubble){
	this(hg,s,t);
	this.firstBubble = true;
    }
    
    //find all ST path in the bubble
    //and remove unsupported paths
    public ArrayList<Path> decompose(Node s, Node t){
	System.err.print("Bubble decomposing...\t");
	this.paths =  this.g.findAllSTPath(s, t);
	System.err.print("Found (" + this.paths.size() + ") possible paths.\n");// + "Removed (");
	
	//resets activePathCoutners in edges
	this.initPathCounters();
	for(Path p : this.paths){
	    //p.printNumActivePaths();
	    p.includePath();
	}
	
	for(Path p : this.paths){
	    p.computeReadSet(g);
	}
	
	//added to keep track of bubble length
	if(this.paths.size() > 0){
	    this.bubbleLengths = new ArrayList<Integer>();
	    this.bubbleLengths.add(new Integer(this.paths.get(0).getPathLength()));
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
    
    public ArrayList<Integer> getBubbleLengths(){
	return this.bubbleLengths;
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

	/* removal of possibly erroneous path */
	/* check for really low weight path compared to other paths in the bubble */
	/* currently using 20% cutoff but we should move to probability model */
	int sumOfReadSetSizeOfSupportedPath = 0;
	int[] readsetSizes = new int[this.paths.size()];

	for(int i=0; i<this.paths.size(); i++){
	    Path p = this.paths.get(i);
	    if(!p.isSupportedPath()){
		readsetSizes[i]=0;
		removalList.add(new Integer(i));
		System.err.print("Removing\tPath" + i + "\t");
		p.printPath();
	    }else{
		readsetSizes[i] = p.getReadSetSize();
		sumOfReadSetSizeOfSupportedPath += readsetSizes[i];
	    }
	}
	
	
	
	for(int i=0; i<readsetSizes.length; i++){
	    if(readsetSizes[i] > 0){
		double ratio = (1.0d * readsetSizes[i]) / ((double) sumOfReadSetSizeOfSupportedPath);
		if( (readsetSizes[i] <= 1 && ratio < 0.2) || 
		    (readsetSizes[i] <= 4 && ratio < 0.1) ||  
		    (readsetSizes[i] > 4 && ratio < 0.05) ){
		    //if(ratio < 0.2){
		    removalList.add(new Integer(i));
		    System.err.print("[Possibly errorneous path] Removing\tPath" + i + "\t");
		    this.paths.get(i).printPath();
		}
	    }
	}
		
	Collections.sort(removalList);

	//removalList is sorted.
	//so we remove from the end so we dont need to worry about index change
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
    
    public MergeStatus mergeBubble(Bubble other, int lastSegregationColumnIndex){
	MergeStatus ms = new MergeStatus(false, lastSegregationColumnIndex);
	int distanceToLastSegregation = other.getStart().get(0) - lastSegregationColumnIndex;
    //public boolean mergeBubble(Bubble other){
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
	    //split first, we will use imputation to connect these spots.
	    if(tpUsed[i] == 0){
		int pathLength = this.getEnd().get(this.getEnd().size()-1) - this.getStart().get(0);
		System.err.println("LOSING TP(" + i + "):\tcurLen:"+pathLength + "\td2ls:" + distanceToLastSegregation);
		if(pathLength >= HLA.READ_LENGTH || distanceToLastSegregation >= (0.5 * HLA.READ_LENGTH)){
		    ms.setSplit(true);
		    System.err.println("CANT PHASE FURTHER. SPLITTING...");
		    return ms;
		}
	    }

	}
	
	
	System.err.println("TOTAL of " + phasedList.size() + "\tphased paths.");
	
	for(int j=0; j<opUsed.length; j++){
	    if(opUsed[j] > 1 && opUsed[j] < phasedList.size()){
		for(int k=0; k<phasedList.size();k++){
		    int[] ijs = phasedList.get(k);
		    if(ijs[1] == j){
			double d = (intersectionSizes.get(k) * 1.0d) / (intersectionSizesSum * 1.0d);
			if( d < 0.1 ){
			    System.err.println("Pruning branch:\td:"+d+"\tTP(" + ijs[0] + ")\tx\tOP(" + ijs[1] + ")");
			    phasedList.remove(k);
			    intersectionSizesSum -= intersectionSizes.get(k);
			    intersectionSizes.remove(k);
			    opUsed[j]--;
			    
			}
		    }
		}
	    }
	}

	int opUsageNum = 0;
	for(int n : opUsed){
	    if(n > 0)
		opUsageNum++;
	}
	
	    
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
	    this.bubbleLengths.addAll(other.getBubbleLengths());
	}
	
	//if we have segregation, meaning we use 2 or more OP(other path)
	if(opUsageNum > 1){
	    ms.setSegregating(true);
	    ms.setLastSegregationColumnIndex(other.getEnd().get(other.getEnd().size() - 1));
	    ms.setSplit(false);
	}
	return ms;
	//return true;
	//}//end else
    }

    public int getNumPaths(){
	return this.paths.size();
    }
    
}
