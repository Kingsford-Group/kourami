import java.io.*;
import java.util.*;
import it.unimi.dsi.fastutil.ints.IntArrayList;

import org.jgrapht.*;
import org.jgrapht.graph.*;

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

    private ArrayList<BubblePathLikelihoodScores> bubbleScores;

    //keeps track of paths through the whole bubble(merged or not-merged)
    private ArrayList<Path> paths;

    private boolean firstBubble;


    public Path getNthPath(int n){
	return this.paths.get(n);
    }
    
    public double getNthBubbleScore(int n, int i, int j){
	if(i<=j)
	    return this.bubbleScores.get(n).getLogScore(i , j);
	else
	    return this.bubbleScores.get(n).getLogScore(j , i);
    }

    public boolean isFirstBubble(){
	return this.firstBubble;
    }

    public int numBubbles(){
	if(this.start.size() == this.bubbleLengths.size())
	    return this.start.size();
	else{
	    System.err.println("---------------->>>>>>>>>>>>>> NOOOOOOOOOOOOOOO!!!!!! ");
	    System.err.println("starts length:\t" + this.start.size());
	    System.err.println("ends length:\t" + this.end.size());
	    System.err.println("bubbleLengths length:\t" + this.bubbleLengths.size());
	    System.err.println("numPaths:\t" + this.paths.size());
	    return -1000;
	}
	
    }

    public void trimPaths(int headerExcess, int tailExcess){
	if(headerExcess > 0 || tailExcess > 0){
	    System.err.println("Trimming by :\t" + headerExcess + "\t" + tailExcess);
	    for(Path p : paths){
		p.trimPath(headerExcess, tailExcess);
	    }
	}
	if(this.bubbleLengths.size() == 1)
	    this.bubbleLengths.set(0, new Integer(this.bubbleLengths.get(0).intValue() - headerExcess - tailExcess));
	else{
	    System.err.println("Something is wrong: Trimming length inconsistent. Exitting");
	    System.exit(-1);
	}
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


    public int printResults(ArrayList<StringBuffer> interBubbleSequences, int startIndex, ArrayList<DNAString> sequences, String hlagenename, int superbubbleNumber, ArrayList<Bubble> bubbles, int bubbleOffset){
	int tmpStartIndex = startIndex;
	
	//for each path (candidate allele)
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
	    /*
	    if(superbubbleNumber == 0 || this.firstBubble){
	    	System.out.print(interBubbleSequences.get(tmpStartIndex).toString());
		curDNA.append(Bubble.stripPadding(interBubbleSequences.get(tmpStartIndex).toString()));
		//output.append(Bubble.stripPadding(interBubbleSequences.get(tmpStartIndex).toString()));
		tmpStartIndex++;
		}*/
	    
	    for(int j=0; j<bubbleSequences.size(); j++){
		if(bubbles.get(j+bubbleOffset).isFirstBubble()){
		    System.err.print("[FB(" + (j+bubbleOffset) +"]");
		    System.out.print(interBubbleSequences.get(tmpStartIndex).toString());
		    curDNA.append(Bubble.stripPadding(interBubbleSequences.get(tmpStartIndex).toString()));
		    tmpStartIndex++;
		}
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
	//for each path in this superbubble
	for(int i=0; i<this.paths.size(); i++){
	    Path p = this.paths.get(i);
	    int pathnum = i;
	    Path curPath = new Path();
	    curPath.setProbability(p.getProbability());
	    curPath.setWeightedIntersectionSum(p.getWeightedIntersectionSum());
	    curPath.setMergedNums(p.getMergedNums());
	    curPath.setMergedTpOpIndicies(p.getMergedTpOpIndicies());
	    curPath.setInterBubbleIntersectionCounts(p.getInterBubbleIntersectionCounts());
	    tmpStartIndex = startIndex;
	    /*if(superbubbleNumber == 0 || this.firstBubble){
		curPath.appendAllEdges(interBubblePaths.get(tmpStartIndex));
		tmpStartIndex++;
		}*/
	    
	    //for each bubble merged in this path
	    int k=0;//k keeps track of index for orderedEdgeList in Path
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

    public int mergePathsInSuperBubbles(ArrayList<TmpPath> interBubblePaths, int startIndex, ArrayList<AllelePath> resultPaths
					, String hlagenename, int superbubbleNumber
					, SimpleDirectedWeightedGraph<Node, CustomWeightedEdge> g, ArrayList<Bubble> bubbles, int bubbleOffset){
	int tmpStartIndex = startIndex;
	
	//for each path in this superbubble: each path here is fractured paths
	for(int i=0; i<this.paths.size(); i++){
	    Path p = this.paths.get(i); 
	    AllelePath curPath = new AllelePath(p.getProbability()
						, p.getWeightedIntersectionSum()
						, p.getMergedNums()
						//, p.getMergedOpIndicies()
						, p);
	
	    ArrayList<ArrayList<CustomWeightedEdge>> bubbleWiseOrderedEdgeLists = p.getBubbleWiseOrderedEdgeList(this.bubbleLengths);
	    //each bubbleWiseOrderedEdgeList is padded by interBubblePaths
	    tmpStartIndex = startIndex;
	    
	    for(int j=0; j<bubbleWiseOrderedEdgeLists.size(); j++){
		// in case of first Bubbles, we need to append interBubble Paths
		if(bubbles.get(j+bubbleOffset).isFirstBubble()){
		    //if(bubbleOffset > 0){//we just mark where disconnectiong in super bubble due to exonic boundaries in AllelePath
		    curPath.setFractureEndIndex();
			//}
		    if(interBubblePaths.get(tmpStartIndex).numEdges() > 0){
			//interBubblePaths.get(tmpStartIndex).toPath(g);
			Path tmpP = interBubblePaths.get(tmpStartIndex).toPath(g);
			curPath.appendAllEdges(tmpP);//interBubblePaths.get(tmpStartIndex).toPath(g));
		    }
		    tmpStartIndex++;
		}
		//then we add bubble edges
		curPath.appendAllEdges(bubbleWiseOrderedEdgeLists.get(j));
		//and then we cap with InterBubble paths
		if(interBubblePaths.get(tmpStartIndex).numEdges() > 0)
		    curPath.appendAllEdges(interBubblePaths.get(tmpStartIndex).toPath(g));
		tmpStartIndex++;
	    }
	    curPath.setSequenceString(g, superbubbleNumber, i);
	    resultPaths.add(curPath);
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
	
	this.bubbleScores = new ArrayList<BubblePathLikelihoodScores>();

	this.decompose(s, t);
	this.removeUnsupported(hg.getGraph(), hg);//this.removeUnsupported();
	
    }

    public Bubble(HLAGraph hg, Node s, Node t, boolean fb){
	this(hg,s,t);
	this.firstBubble = fb;
    }
    
    public Bubble(HLAGraph hg, Node s, Node t, boolean fb, int headerExcessLen, int tailExcessLen){
	this(hg, s, t, fb);
	this.trimPaths(headerExcessLen, tailExcessLen);
    }
    
    //find all ST path in the bubble
    //and remove unsupported paths
    public ArrayList<Path> decompose(Node s, Node t){
	int curBubbleSize = t.getColIndex() - s.getColIndex() + 1;
	//if(curBubbleSize < 20)
	System.err.print("Bubble decomposing...[bubbleSize:" + (t.getColIndex() - s.getColIndex() + 1) +"]\t");
	if(curBubbleSize < 10)
	    this.paths =  this.g.findAllSTPath(s, t);
	else
	    this.paths = this.g.findAllSTPathPruning(s, t);
	
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
    
    public ArrayList<BubblePathLikelihoodScores> getBubbleScores(){
	return this.bubbleScores;
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
    
    public int removeUnsupported(SimpleDirectedWeightedGraph<Node, CustomWeightedEdge> g, HLAGraph hg){
	System.err.println("[Bubble] unsupported path removal...");
	//HashSet<Integer> readHash;
	//CustomHashMap readHash;
	//ArrayList<Integer> removalList = new ArrayList<Integer>();
	IntArrayList removalList = new IntArrayList();
	/* removal of possibly erroneous path */
	/* check for really low weight path compared to other paths in the bubble */
	/* currently using 20% cutoff but we should move to probability model */
	int sumOfReadSetSizeOfSupportedPath = 0;
	int[] readsetSizes = new int[this.paths.size()];
	BubblePathLikelihoodScores scores = null;

	for(int i=0; i<this.paths.size(); i++){
	    Path p = this.paths.get(i);
	    if(!p.isSupportedPath()){
		readsetSizes[i]=0;
		removalList.add(i);//new Integer(i));
		System.err.print("Removing\tPath" + i + "\t");
		p.printPath();
	    }else{
		readsetSizes[i] = p.getReadSetSize();
		sumOfReadSetSizeOfSupportedPath += readsetSizes[i];
	    }
	}

	Collections.sort(removalList);

	/* FIRST we REMOVE paths with no phasing reads */
	//removalList is sorted.
	//so we remove from the end so we dont need to worry about index change
	for(int i=removalList.size() - 1; i >= 0; i--){
	    this.paths.get(removalList.getInt(i)).excludePath();
	    this.paths.remove(removalList.getInt(i));
	}
	System.err.println("Removed (" + removalList.size() + ") paths and\t(" + this.paths.size() + ") left.");
	/* end of phasing reads removal */
	
	removalList = new IntArrayList();

	/* update the readsetSize by copying reassetsize value larger than 0 in order */
	int[] readsetSizeTmp = new int[this.paths.size()];
	int index = 0;
	for(int rs:readsetSizes){
	    if(rs>0){
		readsetSizeTmp[index] = rs;
		index++;
	    }
	}
	readsetSizes = readsetSizeTmp;

	if(HLA.READ_LENGTH < 200){
	    /* First we use a simple heuristic based (parameters) removal */
	    /* Focusing on removing small-weights */
	    for(int i=0; i<readsetSizes.length; i++){
		if(readsetSizes[i] > 0){
		    double ratio = (1.0d * readsetSizes[i]) / ((double) sumOfReadSetSizeOfSupportedPath);
		    if( (readsetSizes[i] <= 1 && ratio < 0.2) || 
			(readsetSizes[i] <= 2 && ratio < 0.134) ||
			(readsetSizes[i] <= 4 && ratio < 0.1) ||  
			(readsetSizes[i] > 4 && ratio < 0.05) ){
			//if(ratio < 0.2){
			removalList.add(i);//new Integer(i));
			System.err.print("[Possibly errorneous path] Removing\tPath" + i + "(|RS|=" +readsetSizes[i] + ",ratio="+ ratio + ",sumReadsetSizes="+ sumOfReadSetSizeOfSupportedPath + ")\t");
			this.paths.get(i).printPath();
		    }
		}
	    }
	}
	
	/* Liklihood calculation retain function  --> geared towards retaining the best */
	if(this.paths.size() >= 2){
	    System.err.println("TEST RUNNING MaxLikelihoodCalcFor BubblePaths");
	    //double[] scoreAndIndices = this.takeMaximumLikeliPair(g);/* [Homo 0-2, Hetero 3-5] 0:homoScore 1:homoIndex1 2:homoIndex2 3:heteroScore 4:heteroIndex1 5: heteroIndex2*/
	    scores = this.takeMaximumLikeliPair(g);
	    double homoScore  = scores.getMaxHomoScore();//scoreAndIndices[0];
	    int homoIndex1 = scores.getMaxHomoGenotypeIndex();//(int) scoreAndIndices[1];
	    int homoIndex2 = homoIndex1;//(int) scoreAndIndices[2];
	    
	    double heteroScore = scores.getMaxHeteroScore();//scoreAndIndices[3];
	    int heteroIndex1 = scores.getMaxHeteroGenotypeIndex1();//(int) scoreAndIndices[4];
	    int heteroIndex2 = scores.getMaxHeteroGenotypeIndex2();//(int) scoreAndIndices[5];
	    int doubleCount1 = scores.getDoubleCountH1();//(int) scoreAndIndices[6];
	    int doubleCount2 = scores.getDoubleCountH1();//(int) scoreAndIndices[7];
	    
	    boolean isAlleleWeak1 = false;
	    boolean isAlleleWeak2 = false;
	    
	    int minCount = doubleCount1;
	    int maxCount = doubleCount2;
	    if(doubleCount1 > doubleCount2){
		minCount = doubleCount2;
		maxCount = doubleCount1;
	    }
	    
	    double frequency = (minCount*1.0d)/ ((minCount+maxCount)*1.0d);
	    if( (minCount <= 2 && frequency < 0.2) ||
		(minCount <= 4 && frequency < 0.134) ||
		(minCount <= 8 && frequency < 0.1) ||
		(minCount  > 8 && frequency < 0.05) ){
		if(doubleCount1 > doubleCount2)
		    isAlleleWeak2 = true;
		else
		    isAlleleWeak1 = true;
	    }
	    
	    double[][] genotypeScores; //i: H1 path index, j: H2 path index

	    System.out.println("Best homozygous genotype is :\t" + homoIndex1 + "/" + homoIndex2 + "\tscore:" + homoScore);
	    System.out.println("Best heterozygous genotype is :\t" + heteroIndex1 + "/" + heteroIndex2 + "\tscore:" + heteroScore);
	    if(isAlleleWeak1)
		System.out.println("H1 count is low : minCount/maxCount" + (minCount/2.0d) + "/" + ((minCount+maxCount)/2.0d));
	    else if(isAlleleWeak2)
		System.out.println("H2 count is low : minCount/maxCount" + (minCount/2.0d) + "/" + ((minCount+maxCount)/2.0d));

	    if(HLA.READ_LENGTH >= 200){
		/*
		if(isAlleleWeak1)
		    removalList.add(heteroIndex1);
		else if(isAlleleWeak2)
		    removalList.add(heteroIndex2);
		
		for(int i=0; i<readsetSizes.length;i++){
		    if( (i != heteroIndex1) && (i != heteroIndex2) ){
			removalList.add(i);
			System.err.println("[Possibly erroneous path] Removing\tPath" + i + "\t");
			this.paths.get(i).printPath();
		    }
		    }*/
		
		double diff = homoScore - heteroScore;
		if(hg.isClassI() && diff > 3){
		    for(int i = 0; i<readsetSizes.length;i++){
			if( i != homoIndex1){
			    removalList.add(i);
			    System.err.println("[Possibly erroneous path] Removing\tPath" + i + "\t");
			    this.paths.get(i).printPath();
			}
		    }
		}else{
		    if(isAlleleWeak1)
			removalList.add(heteroIndex1);
		    else if(isAlleleWeak2)
			removalList.add(heteroIndex2);
		    
		    for(int i=0; i<readsetSizes.length;i++){
			if( (i != heteroIndex1) && (i != heteroIndex2) ){
			    removalList.add(i);
			    System.err.println("[Possibly erroneous path] Removing\tPath" + i + "\t");
			    this.paths.get(i).printPath();
			}
		    }
		    }
	    }else{
		double diff = homoScore - heteroScore;
		System.err.println("diff=homoScore - heteroScore : " + diff);
		if(diff > 0 && diff > 5){ // homoScore > heteroScore
		    /*if(HLA.READ_LENGTH >= 200){	
		      for(int i=0; i<readsetSizes.length;i++){
		      if(i != homoIndex1){
		      removalList.add(i);
		      System.err.println("[Possibly erroneous path] Removing\tPath" + i + "\t");
		      this.paths.get(i).printPath();
		      }
		      }
		      }*/
		    for(int i=0;i<removalList.size();i++){
			if(removalList.getInt(i) == homoIndex1){
			    removalList.removeInt(i);
			    System.err.println("[Retaing path due to liklihood calculation] Rescueing\tPath" + i + "\t");
			}
		    }
		    
		}else{
		    for(int i=0;i<removalList.size();i++){
			if(removalList.getInt(i) == heteroIndex1 || removalList.getInt(i) == heteroIndex2){
			    removalList.removeInt(i);
			    System.err.println("[Retaing path due to liklihood calculation] Rescueing\tPath" + i + "\t");
			}
		    }
		}
	    }	    

	}else{
	    System.err.println("NOT Running likelihood calc because there are less than 2 paths remaining.");
	}


	Collections.sort(removalList);

	//removalList is sorted.
	//so we remove from the end so we dont need to worry about index change
	for(int i=removalList.size() - 1; i >= 0; i--){
	    this.paths.get(removalList.getInt(i)).excludePath();
	    this.paths.remove(removalList.getInt(i));
	}
	System.err.println("Removed (" + removalList.size() + ") paths and\t(" + this.paths.size() + ") left.");
	
	scores.applyRemoval(removalList);
	this.bubbleScores.add(scores);
	return removalList.size();
    }

    
    public PathBaseErrorProb[]  getDataMatrixForLikelihoodCalculation(SimpleDirectedWeightedGraph<Node, CustomWeightedEdge> g){
	PathBaseErrorProb[] pathWisePathBaseErrorProbMatrices = new PathBaseErrorProb[this.paths.size()];
	for(int i=0; i<this.paths.size(); i++)
	    pathWisePathBaseErrorProbMatrices[i] = this.paths.get(i).getBaseErrorProbMatrix(g);
	
	return pathWisePathBaseErrorProbMatrices;
    }
    
    
    //Given all possible paths
    //obtain the pair that gives maxium liklihood of observing data.
    //ONLY considers paths that have phasing-read
    public BubblePathLikelihoodScores takeMaximumLikeliPair(SimpleDirectedWeightedGraph<Node, CustomWeightedEdge> g){
	
	double curScore = 0.0d;
	//obtain path-wise PathBaseErrorProbMaxtrix
	PathBaseErrorProb[] pathWiseErrorProbMatrices = this.getDataMatrixForLikelihoodCalculation(g);

	for(int i=0; i< pathWiseErrorProbMatrices.length; i++){
	    System.err.println("path[" + i + "]:\t");
	    char[] pathBases = pathWiseErrorProbMatrices[i].getBases();
	    for(char c : pathBases){
		System.err.print(c);
	    }
	    System.err.println();
	}

	BubblePathLikelihoodScores scores = new BubblePathLikelihoodScores(this.paths.size());
	
	//getting all possible pairs, including self
	for(int i=0;i<this.paths.size();i++){
	    for(int j=i;j<this.paths.size();j++){
		//NEED TO CALL getScoreForSingleRead over All Reads and sumup the logScore.
		char[] pathBases1 = pathWiseErrorProbMatrices[i].getBases();
		char[] pathBases2 = pathWiseErrorProbMatrices[j].getBases();
		/*if( (i==1 && j==8 ) || (i == 0 && j==1) ){
		    System.out.print("Haplotype[" + i + "] : ");
		    for(char c : pathBases1)
			System.out.print(c);
		    
		    System.out.print("\nHaplotype[" + j + "] : ");
		    for(char c : pathBases2)
			System.out.print(c);
		    System.out.println();
		    }*/
		    
		//iterate over all reads
		curScore = 0.0d;
		Val whichH = new Val();
		int doubleCountH1 = 0;
		int doubleCountH2 = 0;
		for(int k=0; k<this.paths.size();k++){
		    //if( (i== 1 && j == 8) || (i==0 && j == 1) )
		    PathBaseErrorProb curPathErrorMatrix = pathWiseErrorProbMatrices[k];
		    for(int l=0; l<curPathErrorMatrix.numReads(); l++){
			
			boolean debug = false;
			/*
			if( (i== 0 && j == 1) || (i==0 && j == 2) ){
			    char[] readBases = curPathErrorMatrix.getBases();
			    System.out.print("read[" + k + "-" + l + "] : ");
			    for(char c : readBases)
				System.out.print(c);
			    //System.out.println();
			    //debug = true;
			    }*/
			
			double readScore = this.getScoreForSingleRead( curPathErrorMatrix.getNthReadErrorProb(l)
								       , curPathErrorMatrix.getBases()
								       , pathBases1
								       , pathBases2
								       , whichH, debug);// ,i, j);
			/*
			double readScore = this.getScoreForSingleReadMax( curPathErrorMatrix.getNthReadErrorProb(l)
								       , curPathErrorMatrix.getBases()
								       , pathBases1
								       , pathBases2
								       , whichH );// ,i, j);
			*/
			curScore += readScore;
			if( (i== 0 && j == 1) || (i==0 && j == 2) ){
			    System.out.println("\tP(Di|G) = " + readScore);
			}
			if(whichH.getWhichH() == 0)
			    doubleCountH1 += 2;
			else if(whichH.getWhichH() == 1)
			    doubleCountH2 += 2;
			else{
			    doubleCountH1++;
			    doubleCountH2++;
			}
			
		    }
		}
		//System.err.println("logP( D | Haplotype[ " + i + ":" + j + " ] ) =\t" + curScore);
		//System.err.println("logP( D | Haplotype[ " + i + ":" + j + " ] ) =\t" + curScore + "\t|H1|x2=" + doubleCountH1 + "\t|H2|x2=" + doubleCountH2);
		scores.updateMax(i, j, curScore, doubleCountH1, doubleCountH2);
	    }
	}
	return scores;
    }
    
    
    //readPhredScores --> contains ordered phredScore of readBases
    //readBases --> contains ordered readBases {A,C,G,T}
    //pathBases1 --> contains ordered bases {A,C,G,T} of the first allele
    //pathBases2 --> contains ordered bases {A,C,G,T} of the second allele
    private double getScoreForSingleRead(double[] errorProbs, char[] readBases, char[] pathBases1, char[] pathBases2, Val whichH, boolean debug){//, int path_i, int path_j){
	
	//,  columnTransition){
	if(debug){
	    System.err.print("Testing:\t" );
	    for(char c : readBases)
		System.err.print(c);
	    System.err.print("\nAgainst H1:\t");
	    for(char c : pathBases1)
		System.err.print(c);
	    System.err.print("\tH2:\t");
	    for(char c : pathBases2)
		System.err.print(c);
	    System.err.println();
	}
	double logProb1 = 0.0d;
	double logProb2 = 0.0d;
	//for each position in read
	try{
	    for(int i=0; i<readBases.length; i++){
		double errorProb = errorProbs[i];//Math.Pow(0.1, readPhredScores[i]/10.0d);
		double matchProb = 1.0d - errorProb;
		//System.err.println("read(" + i + "):\terroProb:" + errorProb + "\tmatchProb:" + matchProb );
		double mismatchProb = errorProb / 3.0d;
		char readBase = readBases[i];
		char pathBase1 = pathBases1[i];
		char pathBase2 = pathBases2[i];
		if(debug){
		    System.err.println("RB:" + readBase  + "\tPB1:" + pathBases1[i] + "\tPB2:" + pathBases2[i]);
		    System.err.println("MatchP:\t" + matchProb  +"\tMismatchP:\t" + mismatchProb);
		}
		if(readBase == 'N' || pathBase1 == 'N')
		    logProb1 += Math.log(matchProb/4.0d);
		else if(readBase == pathBase1){
		    if(debug)
			System.err.println("MatchPB1");
		    logProb1 += Math.log(matchProb);
		}else{
		    if(debug)
			System.err.println("MismatchPB1");
		    logProb1 += Math.log(mismatchProb);
		}
		if(readBase == 'N' || pathBase2 == 'N')
		    logProb2 += Math.log(matchProb/4.0d);
		else if(readBase == pathBase2){
		    if(debug)
			System.err.println("MatchPB2");
		    logProb2 += Math.log(matchProb);
		}else{
		    if(debug)
			System.err.println("MismatchPB2");
		    logProb2 += Math.log(mismatchProb);
		}
		//System.err.println("logProb1: " + logProb1 + "\tlogProb2: " + logProb2);
		//logProb += Math.log(avgProb);
		//logProb += columnTransition --> need to add transitionProb
	    }
	}catch(ArrayIndexOutOfBoundsException e){
	    System.err.println("|eProbs| :" + errorProbs.length);
	    System.err.println("|rBases| :" + readBases.length);
	    System.err.println("|pBase1| :" + pathBases1.length);
	    System.err.println("|pBase2| :" + pathBases2.length);
	    e.printStackTrace();
	    System.exit(-1);
	}

	/* putting assignment count */
	if(logProb1 > logProb2)
	    whichH.set(0);
	else if(logProb1 == logProb2)
	    whichH.set(-1);
	else
	    whichH.set(1);
	/* end of assignment count*/

	double maxLog = (logProb1 > logProb2 ? logProb1 : logProb2);
	/*
	if( (path_i == 0 && path_j ==1) || (path_i ==1) && (path_j==8) )
	System.err.print("logProb1: " + logProb1 + "\tlogProb2: " + logProb2 + "\tsum: " + (logProb1 + logProb2));
	*/
	//System.err.println("\tlogProb1: " + logProb1 + "\tlogProb2: " + logProb2);
	return maxLog + Math.log(Math.exp(logProb1-maxLog) + Math.exp(logProb2-maxLog)) - Math.log(2);
	
	//return logProb;
    }

    //readPhredScores --> contains ordered phredScore of readBases
    //readBases --> contains ordered readBases {A,C,G,T}
    //pathBases1 --> contains ordered bases {A,C,G,T} of the first allele
    //pathBases2 --> contains ordered bases {A,C,G,T} of the second allele
    private double getScoreForSingleReadMax(double[] errorProbs, char[] readBases, char[] pathBases1, char[] pathBases2, Val whichH){
	
	//,  columnTransition){
	double logProb1 = 0.0d;
	double logProb2 = 0.0d;
	//for each position in read
	try{
	    for(int i=0; i<readBases.length; i++){
		double errorProb = errorProbs[i];//Math.Pow(0.1, readPhredScores[i]/10.0d);
		double matchProb = 1.0d - errorProb;
		double mismatchProb = errorProb / 3.0d;
		char readBase = readBases[i];
		char pathBase1 = pathBases1[i];
		char pathBase2 = pathBases2[i];
		//double avgProb = 0.0d;
		
		
		if(readBase == 'N' || pathBase1 == 'N')
		    logProb1 += Math.log(matchProb/4.0d);
		else if(readBase == pathBase1)
		    logProb1 += Math.log(matchProb);
		else
		    logProb1 += Math.log(mismatchProb);
		
		if(readBase == 'N' || pathBase2 == 'N')
		    logProb2 += Math.log(matchProb/4.0d);
		else if(readBase == pathBase2)
		    logProb2 += Math.log(matchProb);
		else
		    logProb2 += Math.log(mismatchProb);
		
		//logProb += Math.log(avgProb);
		//logProb += columnTransition --> need to add transitionProb
	    }
	}catch(ArrayIndexOutOfBoundsException e){
	    System.err.println("|eProbs| :" + errorProbs.length);
	    System.err.println("|rBases| :" + readBases.length);
	    System.err.println("|pBase1| :" + pathBases1.length);
	    System.err.println("|pBase2| :" + pathBases2.length);
	    e.printStackTrace();
	    System.exit(-1);
	}
	System.err.println("\tlogProb1: " + logProb1 + "\tlogProb2: " + logProb2);
	if(logProb1 > logProb2){
	    whichH.set(0);
	    return logProb1;
	}else if(logProb1 == logProb2){
	    whichH.set(-1);
	    return logProb1;
	}else{
	    whichH.set(1);
	    return logProb2;
	}
    }



    public int removeUnsupported2(){
	System.err.println("[Bubble] unsupported path removal...");
	//HashSet<Integer> readHash;
	CustomHashMap readHash;
	//ArrayList<Integer> removalList = new ArrayList<Integer>();
	IntArrayList removalList = new IntArrayList();
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
		    readHash = e.getReadHashSet().clone();//getReadHashSetDeepCopy();
		    inited = true;
		}else{
		    //update with union after checking intersection
		    //if null, we remove this path
		    

		    if(e.unionAfterCheckingIntersection(readHash) == null){
			removalList.add(i);//new Integer(i));
			break;
		    }
		}
		
	    }
	    System.err.println("Removing\tPath" + i+"\t" + bf.toString());
	    //System.err.print("[" + numUnique +  "]");
	}
	
	for(int i=removalList.size() - 1; i >= 0; i--){
	    this.paths.get(removalList.getInt(i)).excludePath();
	    this.paths.remove(removalList.getInt(i));
	    
	}
	System.err.println("Removed (" + removalList.size() + ") paths and\t(" + this.paths.size() + ") left.");
	return removalList.size();
    }

    /* same as removeUnsupported except this only cares about unique edges */
    public int removeUnsupportedUniqueEdgeOnly(){
	//HashSet<Integer> readHash;
	//ArrayList<Integer> removalList = new ArrayList<Integer>();
	CustomHashMap readHash;
	IntArrayList removalList = new IntArrayList();
	
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
			readHash = e.getReadHashSet().clone();//e.getReadHashSetDeepCopy();
			inited = true;
		    }else{
			//update with union after checking intersection
			//if null, we remove this path
			if(e.unionAfterCheckingIntersection(readHash) == null){
			    removalList.add(i);//new Integer(i));
			    break;
			}
		    }
		}
	    }
	    System.err.print("[" + numUnique +  "]");
	}
	for(int i=removalList.size() - 1; i >= 0; i--){
	    this.paths.get(removalList.getInt(i)).excludePath();
	    this.paths.remove(removalList.getInt(i));
	    
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

    //should only be used between two simple bubble
    //returns int[][] this.paths.sizes() x other.getPaths.size() array. 
    public int[][] getIntersectionCount(Bubble other){
	//int intersectionSizesSum = 0; 
	int[][] intersectionCounts = new int[this.paths.size()][other.getPaths().size()];
	for(int i=0;i<this.paths.size();i++){
	    Path tp = this.paths.get(i);
	    for(int j=0; j<other.getPaths().size();j++){
		Path op = other.getPaths().get(j);
		int intersectionSize = tp.isPhasedWith(op);
		if(intersectionSize >= Path.MIN_SUPPORT_PHASING)
		    intersectionCounts[i][j] = intersectionSize;
		
	    }
	}
	return intersectionCounts;
    }
    
    
    public MergeStatus mergeBubble(Bubble other, int lastSegregationColumnIndex, boolean isClassII, Bubble lastMergedBubble){
	boolean isOtherFirstInInterval = false;
	if(other.isFirstBubble()){
	    isOtherFirstInInterval = true;
	}
	MergeStatus ms = new MergeStatus(false, lastSegregationColumnIndex);
	int distanceToLastSegregation = other.getStart().get(0) - lastSegregationColumnIndex;
    //public boolean mergeBubble(Bubble other){
	ArrayList<Path> paths_new = new ArrayList<Path>();
	//boolean[] tpUsed = new boolean[this.paths.size()];
	//boolean[] opUsed = new boolean[other.getPaths().size()];
	
	int[][] interBubbleIntersectionSizes = lastMergedBubble.getIntersectionCount(other);
	
	/* path used counters */
	int[] tpUsed = new int[this.paths.size()];
	int[] opUsed = new int[other.getPaths().size()];

	
	/* print this paths (DEBUGGIN) */
	for(int i=0;i<this.paths.size();i++){
	    Path tp = this.paths.get(i);
	    System.err.print("TP(" + i + ")\t<readNum:" + tp.getReadSetSize() + ">\t");
	    tp.printInfo();
	    tp.printReadSet();
	}
	
	/* print other paths (DEBUGGIN) */
	for(int i=0;i<other.getPaths().size();i++){
	    Path op = other.getPaths().get(i);
	    System.err.print("OP(" + i + ")\t<readNum:" + op.getReadSetSize() + ">\t");
	    op.printInfo();
	    op.printReadSet();
	}
	
	ArrayList<int[]> phasedList = new ArrayList<int[]>();
	ArrayList<Integer> intersectionSizes = new ArrayList<Integer>();
	int intersectionSizesSum = 0;
	int[] intersectionSizesTPSum = new int[this.paths.size()];
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
		    intersectionSizesTPSum[i] += intersectionSize;
		    intersectionSizes.add(new Integer(intersectionSize));
		    int[] tmp = new int[2];
		    tmp[0] = i;
		    tmp[1] = j;
		    phasedList.add(tmp);
		    tpUsed[i]++;
		    opUsed[j]++;
		}
	    }
	    //decision to split is done after all TP X OP checking
	    /*
	    //split first, we will use imputation to connect these spots.
	    if(tpUsed[i] == 0){
		int pathLength = this.getEnd().get(this.getEnd().size()-1) - this.getStart().get(0);
		System.err.println("LOSING TP(" + i + "):\tcurLen:"+pathLength + "\td2ls:" + distanceToLastSegregation);
		if(pathLength >= HLA.READ_LENGTH || distanceToLastSegregation >= (0.5 * HLA.READ_LENGTH)){
		    ms.setSplit(true);
		    System.err.println("CANT PHASE FURTHER. SPLITTING...");
		    return ms;
		}
		}*/
	}
	
	boolean otherSignificantSignal = false;
	for(Integer i : intersectionSizes){
	    if(i.intValue() >= Path.MIN_SIGNIFICANT_INTERSECTION_SIZE_WHEN_PRUNING)
		otherSignificantSignal = true;
	}
	
	/* Check paths that are being lost due to no phasing reads */
	for(int i=0; i<this.paths.size(); i++){
	    //there is no more phasing reads/pairs
	    if(tpUsed[i] == 0){
		int pathLength = this.getEnd().get(this.getEnd().size()-1) - this.getStart().get(0);
		System.err.println("LOSING TP(" + i + "):\tcurLen:"+pathLength + "\td2ls:" + distanceToLastSegregation);
		if(!otherSignificantSignal){
		    if( (pathLength >= HLA.READ_LENGTH && distanceToLastSegregation >= (0.5 * HLA.READ_LENGTH)) ){
			ms.setSplit(true);
			System.err.println("[FIRST CHECK]CANT PHASE FURTHER. SPLITTING...");
			//return ms;
		    }else if( distanceToLastSegregation >= (1.5 * HLA.READ_LENGTH) ){
			ms.setSplit(true);
			System.err.println("[STRONG SIGNAL]CANT PHASE FURTHER. SPLITTING...");
		    }else if(isClassII && pathLength >=200){
			ms.setSplit(true);
			System.err.println("[CLASS II LENGTH]CANT PHASE FURTHER. SPLITTING...");
		    }
		}else{
		    if( distanceToLastSegregation >= (1.5 * HLA.READ_LENGTH) ){
			ms.setSplit(true);
			System.err.println("[STRONG SIGNAL]CANT PHASE FURTHER. SPLITTING...");
		    }else if(isClassII && pathLength >=200){
			ms.setSplit(true);
			System.err.println("[CLASS II LENGTH]CANT PHASE FURTHER. SPLITTING...");
		    }
		}
	    }
	}
	
	if(ms.isSplit()){
	    for(int i=0; i< this.paths.size(); i++){
		System.err.print("TP(" + i + "):\t");
		this.paths.get(i).printReadSet();
	    }
				 
	    //System.err.println("CANT PHASE FURTHER. SPLITTING...");
	    return ms;
	}
	/* END OF check for no-phasing paths*/
	
	System.err.println("TOTAL of " + phasedList.size() + "\tphased paths.");
	
	int origSizeSum = intersectionSizesSum;

	/*
	for(int j=0; j<opUsed.length; j++){
	    if(opUsed[j] > 1 ){//&& opUsed[j] < phasedList.size()){
		for(int k=0; k<phasedList.size();k++){
		    int[] ijs = phasedList.get(k);
		    if(ijs[1] == j){
			int curSize = intersectionSizes.get(k);
			double d = (intersectionSizes.get(k) * 1.0d) / (origSizeSum * 1.0d);
			System.err.println("Checking branch:\td:" + d + "\tcurSize:" + curSize + "\tTP(" + ijs[0] + ")\tx\tOP(" + ijs[1] + ")");
			//if( (curSize <2 && d < 0.1 ) || (curSize < 3 && d<0.06) || (curSize >= 3 && d <0.05) ){
			if( (curSize <3 && d < 0.1 ) || (curSize >= 3 && d <0.05) ){  
			    System.err.println("Pruning branch:\td:"+d+"\tTP(" + ijs[0] + ")\tx\tOP(" + ijs[1] + ")");
			    phasedList.remove(k);
			    intersectionSizesSum -= intersectionSizes.get(k);
			    intersectionSizes.remove(k);
			    k--;//update the index since we took k out.
			    opUsed[j]--;
			    //added 10/04/16
			    tpUsed[ijs[0]]--;
			    if(tpUsed[ijs[0]] == 0){
				int pathLength = this.getEnd().get(this.getEnd().size()-1) - this.getStart().get(0);
				System.err.println("Pruning results in LOSING TP(" + ijs[0] + "):\tcurLen:"+pathLength + "\td2ls:" + distanceToLastSegregation);
				if(pathLength >= HLA.READ_LENGTH && distanceToLastSegregation >= (0.5 * HLA.READ_LENGTH)){
				    ms.setSplit(true);
				    System.err.println("CANT PHASE FURTHER. SPLITTING... (PRUNING-INDUCED)");
				    return ms;
				}else if(isClassII && pathLength >=200){
				    ms.setSplit(true);
				    System.err.println("CANT PHASE FURTHER. SPLITTING... (CLASS II LENGTH)");
				    return ms;
				}
			    }
			    if(tpUsed[ijs[0]] == 1){
				System.err.println("PRUNING RESULTS IN MORE AGRESSIVE INTERSECTION FOR TP( " + ijs[0] + " )");
			    }
			    //end of added 10/04/16
			}
		    }
		}
	    }
	}
	*/

	int[] tpUsedCopy = tpUsed.clone();
	int[] opUsedCopy = opUsed.clone();

	int origPhasedPathNum = phasedList.size();
	
	for(int i=0; i<phasedList.size();i++){
	    int[] ijs = phasedList.get(i);
	    
	    Path tp = this.paths.get(ijs[0]);
	    int tpUsedCtOrig = tpUsedCopy[ijs[0]];
	    int tpUsedCt = tpUsed[ijs[0]];
	    Path op = other.getPaths().get(ijs[1]);
	    int opUsedCtOrig = opUsedCopy[ijs[1]];
	    int opUsedCt = opUsed[ijs[1]];
	    //either TP or OP is being split.
	    if((tpUsedCtOrig > 1 || opUsedCtOrig > 1)
	       && (tpUsedCt > 0) 
	       && (opUsedCt > 0)){
		int curSize = intersectionSizes.get(i);
		double d = (curSize*1.0d) / (origSizeSum * 1.0d);
		double tpWiseRatio = (curSize*1.0d) / (intersectionSizesTPSum[ijs[0]] * 1.0d);
		System.err.println("Checking branch:\td:" + d + "\ttpWiseRatio:" + tpWiseRatio +  "\tcurSize:" + curSize + "\tTP(" + ijs[0] + ")\tx\tOP(" + ijs[1] + ")"); 

		if( ((curSize <3 && d < 0.1 ) || (curSize >= 3 && d <0.05) || (tpWiseRatio < 0.2)) ){
		    // THIS IS FOR NA12855 ERROR: && origPhasedPathNum > ){ 
		    System.err.println("Pruning branch:\td:"+d+"\tTP(" + ijs[0] + ")\tx\tOP(" + ijs[1] + ")");
		    phasedList.remove(i);
		    intersectionSizesSum -= curSize;
		    intersectionSizes.remove(i);
		    i--;
		    tpUsed[ijs[0]]--;
		    opUsed[ijs[1]]--;
		    if(tpUsed[ijs[0]] == 0){
			int pathLength = this.getEnd().get(this.getEnd().size()-1) - this.getStart().get(0);
			System.err.println("Pruning results in LOSING TP(" + ijs[0] + "):\tcurLen:"+pathLength + "\td2ls:" + distanceToLastSegregation);
			if(pathLength >= HLA.READ_LENGTH && distanceToLastSegregation >= (0.5 * HLA.READ_LENGTH)){
			    ms.setSplit(true);
			    System.err.println("CANT PHASE FURTHER. SPLITTING... (PRUNING-INDUCED)");
			    return ms;
			}else if(isClassII && pathLength >=200){
			    ms.setSplit(true);
			    System.err.println("CANT PHASE FURTHER. SPLITTING... (CLASS II LENGTH)");
			    return ms;
			}
		    }
		    if(tpUsed[ijs[0]] == 1){
			System.err.println("PRUNING RESULTS IN MORE AGRESSIVE INTERSECTION FOR TP( " + ijs[0] + " )");
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
	    paths_new.get(paths_new.size() - 1).updateIntersectionSum(intersectionSize, intersectionSizesSum, ijs, interBubbleIntersectionSizes);
	    
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
	    this.bubbleScores.addAll(other.getBubbleScores());
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


    /* returns indicies for paths that are being phased between two super bubbles 
     * int[0] : path index for this superBubble 
     * int[1] : path index for other supperBubble (sbo)
     * int[2] : intersectionSize
     */
    public ArrayList<int[]> getPhasedSuperBubbles(Bubble sbo){
	
	ArrayList<int[]> phasedList = new ArrayList<int[]>();
	
	//stp: super-this-path
	//sop: super-other-path
	for(int i=0; i<this.paths.size();i++){
	    Path stp = this.paths.get(i);
	    for(int j=0; j < sbo.getPaths().size(); j++){
		Path sop = sbo.getPaths().get(j);
		System.err.print("STP(" + i + ")\tX\tSOP("+j+"):\t");
		int intersectionSize = stp.isPhasedWith(sop);
		//modified so that all stp-sop pair are stored even if they were less the MIN_support_phasing
		//if(intersectionSize >= Path.MIN_SUPPORT_PHASING){
		int[] tmp = new int[3];
		tmp[0] = i;
		tmp[1] = j;
		tmp[2] = intersectionSize;
		phasedList.add(tmp);
		//}
	    }
	}
	return phasedList;
    }
    
    public int getNumPaths(){
	return this.paths.size();
    }
    
}
