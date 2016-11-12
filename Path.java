import it.unimi.dsi.fastutil.ints.IntIterator; 

import htsjdk.samtools.util.QualityUtil;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.Iterator;

import org.jgrapht.*;
import org.jgrapht.graph.*;

public class Path{

    private ArrayList<CustomWeightedEdge> orderedEdgeList;

    private ArrayList<StringBuffer> bubbleSequences; //since path can be used as merged paths, this records bubble sequences in order.
    //private ArrayList<Boolean> isStartBubbles;//true only if the bubble is the first one in a typing interval.
    
    /* adding in for calculating correct prob
    private ArrayList<Double> bubbleIntersectionFraction;
    private ArrayList<int[]> tp-opIndex; //ArrayList of size two int[] {{i_1,j_1},{i_2,j_2}, ... } 
    */
    //private HashSet<Integer> readset;
    /* Key:readID, Value: phredScore --> phredScore on path should not be used. ONLY from rHash in CustomWeightedEdge class */
    // should access phred score from edges in orederedEdgeList.
    private CustomHashMap readset;
    
    public static final int MIN_SUPPORT_BUBBLE = 1;
    
    public static final int MIN_SUPPORT_PHASING = 1;

    public static final int MIN_SIGNIFICANT_INTERSECTION_SIZE_WHEN_PRUNING = 15;

    private double probability;

    private double weightedIntersectionSum;
    
    private int mergedNums;

    private ArrayList<int[]> mergedTpOpIndicies; //keep track of TP-OP index in bubble merging process. Size equals to #of merging done for this path(#bubbles-1)
    private ArrayList<int[][]> interBubbleIntersectionCounts;  //keeps track of counts for all possible tp-op crossing at each merging.
    
    
    public ArrayList<int[][]> getInterBubbleIntersectionCounts(){
	return this.interBubbleIntersectionCounts;
    }
    
    public ArrayList<int[]> getMergedTpOpIndicies(){
	return this.mergedTpOpIndicies;
    }
    
    public void setInterBubbleIntersectionCounts(ArrayList<int[][]> ibic){
	this.interBubbleIntersectionCounts = (ArrayList<int[][]>)ibic.clone();
    }
    
    public void setMergedTpOpIndicies(ArrayList<int[]> moi){
	this.mergedTpOpIndicies = (ArrayList<int[]>)moi.clone();
    }
    
    //this returns TP index that can be used to get corrrect pairing from interBubbleIntersectionCounts.
    public int getLastMergedPathIndex(){
	return this.mergedTpOpIndicies.get(this.mergedTpOpIndicies.size()-1)[1];
    }

    /* THIS IS JUST BASED ON inter-bubble INTERSECTION fraction */
    // NEED TO ADD BUBBLE SPECIFIC GENOTYPE LIKELIHOOD.
    public double[] getJointProbability(Path other, Bubble superBubble){

	double tLogProb = 0.0d;
	double oLogProb = 0.0d;
	double bubblePathLogProb = 0.0d;
	
	//this is the path index from the very first bubble of a superBubble
	int tPreOpIndex = this.mergedTpOpIndicies.get(0)[0];
	int oPreOpIndex = other.getMergedTpOpIndicies().get(0)[0];
	bubblePathLogProb += superBubble.getNthBubbleScore(0, tPreOpIndex, oPreOpIndex);
	
	//for each merging process
	for(int i=0; i<this.mergedTpOpIndicies.size();i++){
	    int tCurOpIndex = this.mergedTpOpIndicies.get(i)[1];
	    int oCurOpIndex = other.getMergedTpOpIndicies().get(i)[1];
	    bubblePathLogProb += superBubble.getNthBubbleScore(i+1, tCurOpIndex, oCurOpIndex);
	    int tSum =  this.sumNthIntersectionCounts(i);
	    int oSum = other.sumNthIntersectionCounts(i);
	    double tFraction = (this.interBubbleIntersectionCounts.get(i)[tPreOpIndex][tCurOpIndex] * 1.0d) / (tSum*1.0d);
	    double oFraction = (other.getInterBubbleIntersectionCounts().get(i)[oPreOpIndex][oCurOpIndex] * 1.0d) / (tSum*1.0d);
	    if(tPreOpIndex == oPreOpIndex 
	       && tCurOpIndex == oCurOpIndex){
		tFraction = tFraction / 2.0d;
		oFraction = oFraction / 2.0d;
	    }
	    tLogProb += Math.log(tFraction);
	    oLogProb += Math.log(oFraction);
	    
	    tPreOpIndex = tCurOpIndex;
	    oPreOpIndex = oCurOpIndex;
	}

	double[] scores = new double[3];//0: intersectionScore for this path, 1: intersectionScore for other path, 2: 
	scores[0] = tLogProb;
	scores[1] = oLogProb;
	scores[2] = bubblePathLogProb;
	return scores;
    }

    public int sumNthIntersectionCounts(int n){
	int[][] tmp = interBubbleIntersectionCounts.get(n);
	int sum = 0;
	for(int[] a : tmp)
	    for(int i : a)
		sum +=i;
	return sum;
    }
    
    /*
    public double[] getJointProbability(Path other){
	double adjustmentLogProb = 0.0d;
	
	int tPre = this.mergedOpIndicies.get(0).intValue();
	int oPre = other.getMergedOpIndicies().get(0).intValue();
	for(int i=1; i<this.mergedOpIndicies.size(); i++){
	    int tCur = this.mergedOpIndicies.get(i).intValue();
	    int oCur = other.mergedOpIndicies.get(i).intValue();
	    if(tPre == oPre && tCur == oCur)
		adjustmentLogProb += Math.log(0.5);
	    
	    tPre = tCur;
	    oPre = oCur;
	}
	
	double[] results = new double[2];
	results[0] = this.probability + adjustmentLogProb;
	results[1] = other.getProbability() + adjustmentLogProb;
	return results;
    }
    
    */
    /*
    public double[] getJointProbability(Path other){
	double adjustmnetLogProb = 0.0d;
	for(int i=0; i < this.mergedTpOpIndicies.size(); i++){
	    int[] tSelection = this.mergedTpOpIndicies.get(i);
	    int[] oSelection = other.getMergedTpOpIndicies().get(i);
	    if(tSelection[0] == oSelection[0] && tSelection[1] == oSelection[1])
		adjustmentLogProb += Math.log(0.5);
	}
	
	double[] results = new double[2];
	results[0] = this.probability + adjustmentLogProb;
	results[1] = other.getProbablity() + adjustmnetLogProb;
	return results;
    }
    */



    //should only be called when path generation via findAllSTPath
    public void trimPath(int headerExcess, int tailExcess){
	if(orderedEdgeList.size() <= (headerExcess + tailExcess)){
	    System.err.println("SERIOUSLY WRONG!!! in trimming header or tail bubble");
	}else{
	    for(int i=0;i<tailExcess;i++)
		this.orderedEdgeList.remove(this.orderedEdgeList.size()-1);
	    
	    for(int j=0; j<headerExcess;j++)
		this.orderedEdgeList.remove(0);
	    
	}
    }


    //breaks the orderedEdgeList as list of each bubble's path.
    public ArrayList<ArrayList<CustomWeightedEdge>> getBubbleWiseOrderedEdgeList(ArrayList<Integer> bubbleLengths){
	ArrayList<ArrayList<CustomWeightedEdge>> bubbleWiseOrderedEdgeList = new ArrayList<ArrayList<CustomWeightedEdge>>();
	int k=0;
	for(int i = 0; i<bubbleLengths.size(); i++){
	    int curlen = bubbleLengths.get(i).intValue();
	    int limit = k + curlen;
	    ArrayList<CustomWeightedEdge> singleBubbleEdgeList = new ArrayList<CustomWeightedEdge>();
	    for(;k<limit;k++){
		try{
		    singleBubbleEdgeList.add(this.orderedEdgeList.get(k));
		}catch(IndexOutOfBoundsException e){
		    System.err.println("curlen=" + curlen + "\tlimit(k+curlen)=" + limit);
		    e.printStackTrace();
		    System.exit(-1);
		}
	    }
	    bubbleWiseOrderedEdgeList.add(singleBubbleEdgeList);
	}
	return bubbleWiseOrderedEdgeList;
    }
    
    //only used during merging interbubbles and superbubbles
    public boolean insertFirstInterBubbleNode(Node n, SimpleDirectedWeightedGraph<Node, CustomWeightedEdge> g){
	//check if edge exists in the graph
	CustomWeightedEdge e = g.getEdge(n, g.getEdgeSource(this.orderedEdgeList.get(0)));
	if(e != null){
	    this.orderedEdgeList.add(0, e);
	    return true;
	}
	return false;
    }


    /* This should ONLY be invoked once to trim down the excess. */
    public void trimExcess(int headerExcess, int tailExcess){
	if(headerExcess > 0){
	    for(int i=0; i<headerExcess; i++)
		this.orderedEdgeList.remove(0);
	    //this.orderedEdgeList.removeRange(0, headerExcess);
	}
	if(tailExcess > 0){
	    for(int i=0; i<tailExcess;i++)
		this.orderedEdgeList.remove(this.orderedEdgeList.size()-1);
	    //this.orderedEdgeList.removeRange(this.orderedEdgeList.size()-tailExcess,this.orderedEdgeList.size());
	}
    }

    public CustomWeightedEdge getNthEdge(int n){
	if(n<this.orderedEdgeList.size())
	    return this.orderedEdgeList.get(n);
	return null;
    }
    
    public void updateIntersectionSum(int intersectionSize, int intersectionSum, int[] mergedTpOpIndex, int[][] interbubbleIntersectionCounts){
	double fraction = (intersectionSize*1.0d)/(intersectionSum*1.0d);
	double logFraction = Math.log(fraction);
	this.probability += logFraction;
	this.weightedIntersectionSum += fraction;
	this.mergedTpOpIndicies.add(mergedTpOpIndex);
	this.interBubbleIntersectionCounts.add(interbubbleIntersectionCounts);
	mergedNums++;
    }

    public PathBaseErrorProb getBaseErrorProbMatrix(SimpleDirectedWeightedGraph<Node, CustomWeightedEdge> g){
	PathBaseErrorProb eProbMatrix = new PathBaseErrorProb(this.readset.size(), this.getPathLength());
	//for each position (edge)
	for(int j=0; j<this.orderedEdgeList.size(); j++){
	    CustomWeightedEdge cur = this.orderedEdgeList.get(j);
	    CustomHashMap curEdgeReadSet = cur.getReadHashSet();
	    char curChar = g.getEdgeTarget(cur).getBase();
	    eProbMatrix.addPathBases(curChar, j);
	    /* for each read */
	    int i = 0;
	    IntIterator itr = this.readset.keySet().iterator();
	    while(itr.hasNext()){
		int curReadID = itr.nextInt();
		int qual = 15;//set it as 15 for unknown 
		if(this.readset.containsKey(curReadID))
		    qual = this.readset.get(curReadID);
		else
		    System.err.println("NO READ COVERAGE");
		//int qual = this.readset.get(curReadID);
		//if(qual == -1)
		//  qual = 15;
		double errorProb = QualityUtil.getErrorProbabilityFromPhredScore(qual);
		eProbMatrix.add(errorProb, i, j);
		i++;
	    }
	}
	
	return eProbMatrix;
    }

    public void printPath(SimpleDirectedWeightedGraph<Node, CustomWeightedEdge> g, int n){//, int headerExcessLen, int tailExcessLen){
	String sequence = this.toString(g, n);//, headerExcessLen, tailExcessLen);
	System.err.println(">candidate_" + n + "\n"+ sequence);
    }

    public String toString(SimpleDirectedWeightedGraph<Node, CustomWeightedEdge> g, int n){//, int headerExcessLen, int tailExcessLen){
	StringBuffer bf = new StringBuffer();
	CustomWeightedEdge pre = null;

	int disconnectCount = 0;
	for(int i=0; i<this.orderedEdgeList.size(); i++){
	    CustomWeightedEdge cur = this.orderedEdgeList.get(i);
	    char curChar;
	    //if(i==0){
	    if(pre == null || !g.getEdgeTarget(pre).equals(g.getEdgeSource(cur))){
		disconnectCount++;
		curChar = g.getEdgeSource(cur).getBase();
		if(curChar != '.')
		    bf.append(curChar);
	    }
	    curChar = g.getEdgeTarget(cur).getBase();
	    if(curChar != '.')
		bf.append(curChar);
	    pre = cur;
	}
	String finalStr = bf.toString();
	/*
	int startIndex = headerExcessLen;
	int endIndex = bf.length() - tailExcessLen;
	String finalStr = bf.substring(startIndex,endIndex);
	*/
	System.err.println(">candidate_" + n + "(" + disconnectCount + ")\n"+ finalStr);//+ bf.toString());
	return finalStr;//bf.toString();
    }


    public void initBubbleSequences(){
	this.bubbleSequences = new ArrayList<StringBuffer>();
	//this.isStartBubbles = new ArrayList<Boolean>();
    }

    //Should only be used when it's NOT a merged bubble
    public int getPathLength(){
	return this.orderedEdgeList.size();
    }
    
    public int getMergedNums(){
	return this.mergedNums;
    }

    public double getProbability(){
	return this.probability;
    }

    public double getWeightedIntersectionSum(){
	return this.weightedIntersectionSum;
    }

    public double getAvgWeightedIntersectionSum(){
	return this.weightedIntersectionSum/mergedNums;
    }

    public void setProbability(double p){
	this.probability = p;
    }

    public void setWeightedIntersectionSum(double s){
	this.weightedIntersectionSum = s;
    }
    
    public void setMergedNums(int n){
	this.mergedNums = n;
    }
    
    public CustomHashMap getReadSet(){
	return this.readset;
    }

    public void printReadSet(){
	this.readset.printReads();
    }

    public int getReadSetSize(){
	return this.readset.size();
    }

    public ArrayList<StringBuffer> getBubbleSequences(){
	return this.bubbleSequences;
    }
    
    /*
    public ArrayList<StringBuffer> getIsStartBubbles(){
	return this.isStartBubbles;
    }
    */

    public void setReadSet(CustomHashMap rs){//HashSet<Integer> rs){
	this.readset = rs;
    }

    public void subtractReadSet(CustomHashMap ors){//HashSet<Integer> ors){
	this.readset.removeAll(ors);
    }

    public void printBubbleSequence(SimpleDirectedWeightedGraph<Node, CustomWeightedEdge> g){
	StringBuffer bf = new StringBuffer();
	for(int i=0;i<this.orderedEdgeList.size()-1;i++){
	    CustomWeightedEdge e = this.orderedEdgeList.get(i);
	    CustomWeightedEdge e2 = this.orderedEdgeList.get(i+1);
	    bf.append(g.getEdgeTarget(e).getBase() + "(" + g.getEdgeWeight(e)+":"+g.getEdgeWeight(e2)+")");
	}
	System.err.println("bubble:\t" + bf.toString());
    }

    public void initBubbleSequence(SimpleDirectedWeightedGraph<Node, CustomWeightedEdge> g){
	if(bubbleSequences.size() == 0){
	    StringBuffer bf = new StringBuffer();
	    for(int i=0;i<this.orderedEdgeList.size()-1;i++){
		CustomWeightedEdge e = this.orderedEdgeList.get(i);
		bf.append(g.getEdgeTarget(e).getBase());
	    }
	    System.err.println("bubble:\t" + bf.toString());
	    bubbleSequences.add(bf);
	}else{
	    System.err.println("Shouldn't be called here.");
	    System.exit(-1);
	}
    }
    
    public void addBubbleSequence(StringBuffer sb){
	this.bubbleSequences.add(sb);
    }

    public Path deepCopy(){
	Path p = new Path();
	for(CustomWeightedEdge e : this.orderedEdgeList){
	    p.appendEdge(e);
	}
	for(StringBuffer sb : this.bubbleSequences){
	    p.addBubbleSequence(sb);
	}
	p.setReadSet(this.readset.clone());//this.getReadSetDeepCopy());
	
	p.setProbability(this.probability);
	p.setWeightedIntersectionSum(this.weightedIntersectionSum);
	p.setMergedNums(this.mergedNums);
	p.setMergedTpOpIndicies(this.mergedTpOpIndicies);
	p.setInterBubbleIntersectionCounts(this.interBubbleIntersectionCounts);
	return p;
    }
    

    public boolean isSupportedPath(){
	if(readset.size() >= Path.MIN_SUPPORT_BUBBLE){
	    return true;
	}
	return false;
    }

    public void computeReadSet(HLAGraph g){
	System.err.println("Verifying:");
	this.printPath(g);
	//HashSet<Integer> tmpset = new HashSet<Integer>();
	//HashSet<Integer> unionUniqueSet = new HashSet<Integer>();
	CustomHashMap tmpset = new CustomHashMap();
	CustomHashMap unionUniqueSet = new CustomHashMap();
	
	//first check if size of intersection is nonzero.
 	for(int i=0; i<this.orderedEdgeList.size(); i++){
	    CustomWeightedEdge e = this.orderedEdgeList.get(i);
	    if(e.isUniqueEdge())
		unionUniqueSet.addAll(e.getReadHashSet());
	    
	    if(i == 0)
		tmpset.addAll(e.getReadHashSet());
	    else{
		if(tmpset.size() > 0)
		    tmpset.intersectionPE(e.getReadHashSet());//we take intersection
		//no need to break in case there are unique edges.
		//if(tmpset.size() == 0)
		//    break;
	    }
	}
	
	//intersection is nonzero, we will add supplemnentary evidences(uniqEdgeReads union)
	if(tmpset.size() >= Path.MIN_SUPPORT_BUBBLE ){
	    System.err.print("InersectionSize\t" + tmpset.size()+ "\tUnionUniqSetSize\t" + unionUniqueSet.size());
	    /*
	    ArrayList<CustomWeightedEdge> nonUniqueEdges = new ArrayList<CustomWeightedEdge>();
	    for(CustomWeightedEdge e : this.orderedEdgeList){
		if(!e.isUniqueEdge()){
		    nonUniqueEdges.add(e);
		    HashSet<Integer> tmpset2 = new HashSet<Integer>();
		    tmpset2.addAll(e.getReadHashSet());
		    tmpset2.retainAll(unionUniqueSet);
		}
		}*/
	    tmpset.addAll(unionUniqueSet);
	    System.err.println("\tTotalSetSize\t" + tmpset.size());
	    for(CustomWeightedEdge e : this.orderedEdgeList){
		e.subtractSet(tmpset);
	    }
	}else{
	    System.err.print("InersectionSize\t" + tmpset.size()+ "\tUnionUniqSetSize\t" + unionUniqueSet.size());
	    tmpset = new CustomHashMap();//new HashSet<Integer>();
	    System.err.println("TotalSetSize\t" + tmpset.size() + "\t----> REMOVED");
	}
	this.readset = tmpset;
    }


    
    //intersection of reads acorss edges
    //also add reads belonging to unique edges
    //this works because the bubble is small.
    /*public void computeReadSet(){
	for(int i=0; i<this.orderedEdgeList.size(); i++){
	    CustomWeightedEdge e = this.orderedEdgeList.get(i);
	    if(i==0){
		this.readset = e.getReadHashSetDeepCopy();
	    }else{
		this.readset.retainAll(e.getReadHashSet());
		if(this.readSet.size() == 0){
		    break;
		}
	    }
	}
	
	if(this.readSet.size() > 0){
	    for(CustomWeightedEdge e : this.orderedEdgeList)
		e.substractSet(this.readset);
	}
    }
    */
    

    public ArrayList<CustomWeightedEdge> getOrderedEdgeList(){
	return this.orderedEdgeList;
    }

    public Path(){
	this.orderedEdgeList = new ArrayList<CustomWeightedEdge>();
	this.readset = new CustomHashMap();//new HashSet<Integer>();
	this.bubbleSequences = new ArrayList<StringBuffer>();
	this.weightedIntersectionSum = 0.0d;
	this.mergedNums = 0;
	this.probability = 0.0d;
	this.mergedTpOpIndicies = new ArrayList<int[]>();
	this.interBubbleIntersectionCounts = new ArrayList<int[][]>();
    }

    public Path(double p, double wis, int mn){
	this();
	this.weightedIntersectionSum = wis;
	this.mergedNums = mn;
	this.probability = p;
    }


    public void appendEdge(CustomWeightedEdge e){
	this.orderedEdgeList.add(e);
    }
    
    public void appendAllEdges(Path other){
	this.orderedEdgeList.addAll(other.getOrderedEdgeList());
    }

    public void appendAllEdges(ArrayList<CustomWeightedEdge> anotherOrderedEdgeList){
	this.orderedEdgeList.addAll(anotherOrderedEdgeList);
    }


    public Path combinePaths(Path other){
	Path np = this.deepCopy();
	np.appendAllEdges(other);
	np.setReadSet(new CustomHashMap());//new HashSet<Integer>());
	np.initBubbleSequences();
	np.setWeightedIntersectionSum(np.getWeightedIntersectionSum() + other.getWeightedIntersectionSum());
	np.setMergedNums(np.getMergedNums() + other.getMergedNums());
	np.setProbability(np.getProbability() + other.getProbability());
	return np;
    }

    public Path(CustomWeightedEdge e){
	this();
	this.appendEdge(e);
    }

    public Node getLastVertex(SimpleDirectedWeightedGraph<Node, CustomWeightedEdge> g){
	if(this.orderedEdgeList == null || this.orderedEdgeList.size() == 0)
	    return null;
	else{
	    return g.getEdgeTarget(this.orderedEdgeList.get(this.orderedEdgeList.size() - 1));
	}
    }
    
    // M : M 
    // tp is used multiple times and op is used multiple times
    public Path mergePathManytoMany(Path other){
	Path np = this.mergePaths(other);
	np.getIntersectionPE(other); // replaced with paired-end aware intersection
	return np;
    }

    // 1 : M
    //tp is used once but op is used multiple times
    public Path mergePath1toMany(Path other){
	Path np = this.mergePaths(other);
	//updated to use intersectionPE. instead 
	this.readset.addPEReads(other.getReadSet());
	//np.setReadSet(this.readset.clone().union(this.readset.clone().intersectionPE(other.getReadSet())));
	return np;
    }

    // M : 1
    //tp is used multiple times but op is used once.
    public Path mergePathManyto1(Path other){
	Path np = this.mergePaths(other);
	CustomHashMap tmp = other.getReadSet().clone();
	tmp.addPEReads(this.readset);
	np.setReadSet(tmp);
       
	//np.setReadSet(other.getReadSet().clone());//other.getReadSetDeepCopy());
	//np.setReadSet(other.getReadSet().clone().union(other.getReadSet().clone().intersectionPE(this.getReadSet())));
	return np;
    }

    // 1 : 1
    //1 to 1 connection tp and op are not split into other flows
    public Path mergePathsUnique(Path other){
	Path np = this.mergePaths(other);
	np.unionReadSets(other);
	return np;
    }
    
    private void unionReadSets(Path other){
	this.readset.addAll(other.getReadSet());
    }
    
    /* this simply merged edges */
    /* readSet is handled separately */
    public Path mergePaths(Path other){
	Path p = this.deepCopy();
	p.getOrderedEdgeList().addAll(other.getOrderedEdgeList());
	p.getBubbleSequences().addAll(other.getBubbleSequences());
	return p;
    }

    /* phasing based on unique edges */
    //checks for phasing support information between two paths
    //Two paths are phased if intersection of reads passing through unique edges.
    //@returns true if intersecting sets are NOT empty
    //@returns false otherwise.
    public boolean isPhasedWithOLD(Path other){
	//this set
	//HashSet<Integer> ts = this.getUnionOfUniqueEdgesReadSet();
	CustomHashMap ts = this.getUnionOfUniqueEdgesReadSet();
	//other set
	//HashSet<Integer> os = other.getUnionOfUniqueEdgesReadSet();
	CustomHashMap os = other.getUnionOfUniqueEdgesReadSet();
	System.out.println("TS:\t");
	ts.printKeys();//Path.printHashSet(ts);
	System.out.println("OS:\t");
	os.printKeys();//Path.printHashSet(os);
	
	ts.intersectionPE(os);
	System.out.print("\t");
	ts.printKeys();//Path.printHashSet(ts);
	
	if(ts.size() >= Path.MIN_SUPPORT_PHASING)
	    return true;
	return false;
    }
    
    public void getIntersection(Path other){
	this.readset.retainAll(other.getReadSet());
    }
    
    public void getIntersectionPE(Path other){
	this.readset.intersectionPE(other.getReadSet());
    }

    //NO LONGER USED. SHOULD USE clone() in CustomHashMap class.
    /*
    private HashSet<Integer> getReadSetDeepCopy(){
    	HashSet<Integer> tmp = new HashSet<Integer>();
	
	Iterator<Integer> itr = this.readset.iterator();
	while(itr.hasNext()){
	    tmp.add(itr.next());
	}
	return tmp;
	}*/

    /* phasing based on read set */
    public int isPhasedWith(Path other){
	//HashSet<Integer> copyset = this.getReadSetDeepCopy();
	CustomHashMap copyset = this.readset.clone();
	//copyset.retainAll(other.getReadSet());
	copyset.intersectionPE(other.getReadSet());// special intersection for paired-end 
	if(copyset.size() >= Path.MIN_SUPPORT_PHASING){
	    System.err.println("PHASED[intersectionSize:" + copyset.size() + "]");
	    //return true;
	}else{
	    System.err.println("NOT PHASED[intersectionSize:" + copyset.size() + "]");
	    if(copyset.size() > 0){
		this.subtractReadSet(copyset);
		other.subtractReadSet(copyset);
	    }
	    //return false;
	}
	return copyset.size();
    }

    public String getNumUniqueEdges(){
	int count = 0;
	StringBuffer bf = new StringBuffer();
	for(CustomWeightedEdge e : this.orderedEdgeList){
	    if(e.isUniqueEdge()){
		bf.append("{"+e.getEdgeId()+"}");
		count++;
	    }
	}
	bf.append("(" + count + ")");
	return bf.toString();
	//return count;
    }
    
    public void printInfo(){
	System.err.print("NumEdges:" + this.orderedEdgeList.size());
	System.err.print("\tNumUniqueEdges:" + this.getNumUniqueEdges() + "\n");
    }

    public String toSimplePathString(HLAGraph g){
	StringBuffer bf = new StringBuffer();
	for(CustomWeightedEdge e : this.orderedEdgeList)
	    bf.append(g.getGraph().getEdgeTarget(e).getBase());
	return bf.toString();
    }

    public void printPath(){
	System.err.print("NumEdges:" + this.orderedEdgeList.size() + "\t");
	for(CustomWeightedEdge e : this.orderedEdgeList){
	    System.err.print("{"+e.getEdgeId()+"}");
	}
	System.err.println();
    }

    public void printPath(HLAGraph g){
	System.err.print("NumEdges:" + this.orderedEdgeList.size() + "\t");

	for(CustomWeightedEdge e : this.orderedEdgeList){
	    System.err.print("{"+e.getEdgeId()+"}" + g.getGraph().getEdgeTarget(e).getBase());
	}
	System.err.println();
    }

    
    //if there is no uniqueEdges, return null
    //else it returns union HashSet<Integer> of reads over all unique edges.
    //size 0 if there is reads covering unique edge
    //    public HashSet<Integer> getUnionOfUniqueEdgesReadSet(){
    public CustomHashMap getUnionOfUniqueEdgesReadSet(){
	System.err.println("UnionOfUniqueEdges");
	CustomHashMap s = new CustomHashMap();
	boolean atLeast1UniqueEdge = false;
	for(CustomWeightedEdge e : this.orderedEdgeList){
	    System.err.print("|" + e.getNumActivePath() + "|");
	    if(e.isUniqueEdge()){
		System.err.println("U:"+ e.getEdgeId());
		atLeast1UniqueEdge  = true;
		s.addAll(e.getReadHashSet());
	    }else
		System.err.println("R:"+ e.getEdgeId());
	}
	if(!atLeast1UniqueEdge)
	    return null;
	return s;
    }
    
    public void initPathCounter(){
	for(CustomWeightedEdge e : this.orderedEdgeList){
	    e.initNumActivePath(); //sets it to zero 
	    //e.includeEdge(); // increment numActivePath
	}
    }

    public void includePath(){
	for(CustomWeightedEdge e : this.orderedEdgeList){
	    e.includeEdge();
	}
    }
    
    public void includePathNTimes(int n){
	for(CustomWeightedEdge e : this.orderedEdgeList){
	    e.includeEdgeNTimes(n);
	}
    }

    public void excludePath(){
	for(CustomWeightedEdge e : this.orderedEdgeList){
	    e.excludeEdge();
	}
    }

    public void printNumActivePaths(){
	for(CustomWeightedEdge e : this.orderedEdgeList){
	    System.err.print(e.getEdgeId() + "\t" + e.getNumActivePath() + "\t|\t");
	}
	System.err.println();
    }
}
