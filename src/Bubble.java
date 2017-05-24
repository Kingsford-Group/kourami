/*
Part of Kourami HLA typer/assembler
(c) 2017 by  Heewook Lee, Carl Kingsford, and Carnegie Mellon University.
See LICENSE for licensing.
*/
import java.util.*;
import it.unimi.dsi.fastutil.ints.IntArrayList;
import it.unimi.dsi.fastutil.ints.IntIterator;  

import org.jgrapht.*;
import org.jgrapht.graph.*;

public class Bubble{

    private HLAGraph g;

    private ArrayList<Integer> start;
    private ArrayList<Integer> end;
    
    //s- and t- nodes of each bubbles that are merged
    private ArrayList<Node> sNodes;
    private ArrayList<Node> tNodes;

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
    
    public double getNthBubbleFractionScore(int n, int i, int j){
	if(i<=j)
	    return this.bubbleScores.get(n).getLogFractionScore(i , j);
	else
	    return this.bubbleScores.get(n).getLogFractionScore(j , i);
    }

    public boolean isFirstBubble(){
	return this.firstBubble;
    }

    public int numBubbles(){
	if(this.start.size() == this.bubbleLengths.size())
	    return this.start.size();
	else{
	    if(HLA.DEBUG){
		HLA.log.appendln("---------------->>>>>>>>>>>>>> NOOOOOOOOOOOOOOO!!!!!! ");
		HLA.log.appendln("starts length:\t" + this.start.size());
		HLA.log.appendln("ends length:\t" + this.end.size());
		HLA.log.appendln("bubbleLengths length:\t" + this.bubbleLengths.size());
		HLA.log.appendln("numPaths:\t" + this.paths.size());
	    }
	    return -1000;
	}
	
    }
    
    public void trimPaths(int headerExcess, int tailExcess){
	if(headerExcess > 0 || tailExcess > 0){
	    if(HLA.DEBUG)
		HLA.log.appendln("Trimming by :\t" + headerExcess + "\t" + tailExcess);
	    for(Path p : paths){
		p.trimPath(headerExcess, tailExcess);
	    }
	}
	if(this.bubbleLengths.size() == 1)
	    this.bubbleLengths.set(0, new Integer(this.bubbleLengths.get(0).intValue() - headerExcess - tailExcess));
	else{
	    HLA.log.appendln("Something is wrong: Trimming length inconsistent. Exitting");
	    HLA.log.outToFile();
	    System.exit(-1);
	}
    }

    public void printBubbleSequenceSizes(){
	for(Path p : this.paths)
	    HLA.log.append(p.getBubbleSequences().size() + "\t");
	HLA.log.appendln();
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
	HLA.log.appendln("Printing\t" + this.paths.size() + "\tpossible sequences"  );
	for(int i=0; i<this.paths.size(); i++){
	    Path p = this.paths.get(i);
	    ArrayList<StringBuffer> bubbleSequences = p.getBubbleSequences();
	    HLA.log.appendln(">>>>>>>ITBS\t" + interBubbleSequences.size() + "\tBS\t" + bubbleSequences.size());
	    HLA.log.append(interBubbleSequences.get(0).toString());
	    for(int j=1; j<interBubbleSequences.size(); j++){
		HLA.log.append("<<" + bubbleSequences.get(j-1).toString() + ">>");
		HLA.log.append(interBubbleSequences.get(j).toString());
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
	    HLA.log.appendln("IntersectionScore:\t" + p.getAvgWeightedIntersectionSum() + "\t" + p.getProbability());
	    ArrayList<StringBuffer> bubbleSequences = p.getBubbleSequences();
	    //each bubbleSequence is padded by interBubbleSequences
	    //so we print the first interBubbleSequence.
	    tmpStartIndex = startIndex;
	    if(superbubbleNumber == 0 || this.firstBubble){
		HLA.log.append(interBubbleSequences.get(tmpStartIndex).toString());
		curDNA.append(Bubble.stripPadding(interBubbleSequences.get(tmpStartIndex).toString()));
		//output.append(Bubble.stripPadding(interBubbleSequences.get(tmpStartIndex).toString()));
		tmpStartIndex++;
	    }
	    
	    
	    for(int j=0; j<bubbleSequences.size(); j++){
		HLA.log.append(" <" + bubbleSequences.get(j) + "> ");//prints the bubble
		curDNA.append(Bubble.stripPadding(bubbleSequences.get(j).toString()));
		//output.append(Bubble.stripPadding(bubbleSequences.get(j).toString()));
		//if path ends with bubble, we dont need to print the last interBubbleSequence
		if(tmpStartIndex < interBubbleSequences.size()){
		    HLA.log.append(interBubbleSequences.get(tmpStartIndex).toString()); //prints the interBubble
		    curDNA.append(Bubble.stripPadding(interBubbleSequences.get(tmpStartIndex).toString()));
		    //output.append(Bubble.stripPadding(interBubbleSequences.get(tmpStartIndex).toString()));
		}
		
		tmpStartIndex++;
	    }
	    HLA.log.appendln();
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

	    HLA.log.appendln("IntersectionScore:\t" + p.getAvgWeightedIntersectionSum() + "\t" + p.getProbability());
	    ArrayList<StringBuffer> bubbleSequences = p.getBubbleSequences();
	    //each bubbleSequence is padded by interBubbleSequences
	    //so we print the first interBubbleSequence.
	    tmpStartIndex = startIndex;
	    
	    for(int j=0; j<bubbleSequences.size(); j++){
		if(bubbles.get(j+bubbleOffset).isFirstBubble()){
		    HLA.log.append("[FB(" + (j+bubbleOffset) +"]");
		    HLA.log.append(interBubbleSequences.get(tmpStartIndex).toString());
		    curDNA.append(Bubble.stripPadding(interBubbleSequences.get(tmpStartIndex).toString()));
		    tmpStartIndex++;
		}
		HLA.log.append(" <" + bubbleSequences.get(j) + "> ");//prints the bubble
		curDNA.append(Bubble.stripPadding(bubbleSequences.get(j).toString()));

		//if path ends with bubble, we dont need to print the last interBubbleSequence
		if(tmpStartIndex < interBubbleSequences.size()){
		    HLA.log.append(interBubbleSequences.get(tmpStartIndex).toString()); //prints the interBubble
		    curDNA.append(Bubble.stripPadding(interBubbleSequences.get(tmpStartIndex).toString()));
		}
		tmpStartIndex++;
	    }

	    HLA.log.appendln();
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
	    curPath.setInterBubbleIntersectionCumulativeCounts(p.getInterBubbleIntersectionCumulativeCounts());
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

    public Bubble(HLAGraph hg, Node s, Node t){
	this(hg, s, t, null, null);
    }
    

    public Bubble(HLAGraph hg, Node s, Node t, Node[] headerNodes, Node[] tailNodes){
	this(hg, s, t, false, headerNodes, tailNodes);
    }
    
    public Bubble(HLAGraph hg, Node s, Node t, boolean fb, Node[] headerNodes, Node[] tailNodes){
	this.firstBubble = fb;
	this.g = hg;
	this.sNodes = new ArrayList<Node>();
	this.tNodes = new ArrayList<Node>();
	this.sNodes.add(s);
	this.tNodes.add(t);
	this.start = new ArrayList<Integer>();
	this.end = new ArrayList<Integer>();
	if(s == null){
	    HLA.log.appendln("[BUBBLE] start node null");
	}
	this.start.add(new Integer(s.getColIndex()));
	this.end.add(new Integer(t.getColIndex()));
	this.paths = new ArrayList<Path>();
	
	this.bubbleScores = new ArrayList<BubblePathLikelihoodScores>();
	
	this.decompose(s, t, headerNodes, tailNodes);
	this.removeUnsupported(hg.getGraph(), hg, s, t);//this.removeUnsupported();
    }
    
    public Bubble(HLAGraph hg, Node s, Node t, boolean fb, int headerExcessLen, int tailExcessLen, Node[] headerNodes, Node[] tailNodes){
	this(hg, s, t, fb, headerNodes, tailNodes);
	//this.trimPaths(headerExcessLen, tailExcessLen);
    }
  
    private void updateStart(Node newS){
	if(this.sNodes.size() == 1){
	    this.sNodes.set(0, newS);
	    this.start.set(0, newS.getColIndex());
	}
    }
    
    private void updateEnd(Node newE){
	if(this.tNodes.size() == 1){
	    this.tNodes.set(0, newE);
	    this.end.set(0, newE.getColIndex());
	}
    }

  
    //find all ST path in the bubble
    //and remove unsupported paths
    public ArrayList<Path> decompose(Node s, Node t, Node[] headerNodes, Node[] tailNodes){
	int curBubbleSize = t.getColIndex() - s.getColIndex() + 1;
	//if(curBubbleSize < 20)
	if(HLA.DEBUG){
	    if(headerNodes == null && tailNodes == null){
		HLA.log.append("Bubble decomposing...[bubbleSize:" + (t.getColIndex() - s.getColIndex() + 1) +"]\t");
	    }
	    if(headerNodes != null){
		for(Node n: headerNodes)
		    HLA.log.appendln("HN:\t" + n.toString());
	    }
	    if(tailNodes != null){
		for(Node n : tailNodes)
		    HLA.log.appendln("TN:\t" + n.toString());
	    }
	}
	
	if(headerNodes == null && tailNodes == null){
	    if(curBubbleSize < 10)
		this.paths = this.g.findAllSTPath(s, t);
	    else
		this.paths = this.g.findAllSTPathPruning(s, t);
	}else if(headerNodes == null && tailNodes !=null){
	    curBubbleSize = tailNodes[0].getColIndex() - s.getColIndex() + 1;
	    if(HLA.DEBUG){
		HLA.log.append("[T]Bubble decomposing...[bubbleSize:" + curBubbleSize +"]\t");
		HLA.log.appendln(s.toString() + ":");
	    }
	    for(Node tn : tailNodes){
		if(HLA.DEBUG)
		    HLA.log.appendln("\t" + tn.toString());
		if(curBubbleSize < 10)
		    this.paths.addAll(this.g.findAllSTPath(s, tn));
		else
		    this.paths.addAll(this.g.findAllSTPathPruning(s, tn));
	    }
	    updateEnd(tailNodes[0]);
	}else if(headerNodes !=null && tailNodes == null){
	    curBubbleSize = t.getColIndex() - headerNodes[0].getColIndex() + 1;
	    if(HLA.DEBUG)
		HLA.log.append("[H]Bubble decomposing...[bubbleSize:" + curBubbleSize +"]\t");
	    for(Node sn : headerNodes){
		if(curBubbleSize < 10)
		    this.paths.addAll(this.g.findAllSTPath(sn, t));
		else
		    this.paths.addAll(this.g.findAllSTPathPruning(sn, t));
	    }
	    updateStart(headerNodes[0]);
	}else{
	    curBubbleSize = tailNodes[0].getColIndex() - headerNodes[0].getColIndex() + 1;
	    if(HLA.DEBUG)
		HLA.log.append("[HT]Bubble decomposing...[bubbleSize:" + curBubbleSize +"]\t");
	    for(Node sn : headerNodes){
		for(Node tn : tailNodes){
		    if(curBubbleSize > 10)
			this.paths.addAll(this.g.findAllSTPath(sn, tn));
		    else
			this.paths.addAll(this.g.findAllSTPathPruning(sn, tn));
		}
	    }
	}
	
	if(HLA.DEBUG)
	    HLA.log.append("Found (" + this.paths.size() + ") possible paths.\n");// + "Removed (");
	
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
	    HLA.log.append("\tP" + count + "\t");
	    p.printPath();
	    count++;
	}
    }
    
    public int removeUnsupported(SimpleDirectedWeightedGraph<Node, CustomWeightedEdge> g, HLAGraph hg, Node s, Node t){
	if(HLA.DEBUG)
	    HLA.log.appendln("[Bubble] unsupported path removal...");
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
	int bubbleSize = t.getColIndex() - s.getColIndex() + 1;
	
	for(int i=0; i<this.paths.size(); i++){
	    Path p = this.paths.get(i);
	    int numEdges = p.getOrderedEdgeList().size();
	    /*if(numEdges != (bubbleSize - 1)){
		readsetSizes[i]=0;
		removalList.add(i);//new Integer(i));
		if(HLA.DEBUG){
		    HLA.log.append("Removing(WRONG PATH LENGTH)\tPath" + i + "\t");
		    p.printPath();
		}
	    }
	    else */if(!p.isSupportedPath()){
		readsetSizes[i]=0;
		removalList.add(i);//new Integer(i));
		if(HLA.DEBUG3){
		    HLA.log.append("Removing\tPath" + i + "\t");
		    p.printPath();
		}
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
	if(HLA.DEBUG)
	    HLA.log.appendln("Removed (" + removalList.size() + ") paths and\t(" + this.paths.size() + ") left.");
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
	if(HLA.DEBUG){
	    if(this.isFirstBubble()){
		HLA.log.appendln(">>>>>FIRSTBUBBLE<<<<<<");
	    }
	}
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
			if(HLA.DEBUG3){
			    HLA.log.append("[Possibly errorneous path] Removing\tPath" + i + "(|RS|=" +readsetSizes[i] + ",ratio="+ ratio + ",sumReadsetSizes="+ sumOfReadSetSizeOfSupportedPath + ")\t");
			    this.paths.get(i).printPath();
			}
		    }
		}
	    }
	}else{
	    if(this.isFirstBubble()){
		for(int i=0; i<readsetSizes.length; i++){
		    if(readsetSizes[i] > 0){
			if(HLA.DEBUG3)
			    HLA.log.appendln("===== CHECKING PATH[" + i + "]");
			double ratio = (1.0d * readsetSizes[i]) / ((double) sumOfReadSetSizeOfSupportedPath);
			if(ratio < 0.2){
			    if(HLA.DEBUG3)
				HLA.log.appendln("===== ratio:\t" + ratio + " =======");
			    CustomHashMap tMap = this.paths.get(i).getReadSet();
			    for(int j=0;j<readsetSizes.length;j++){
				if(i!=j){
				    double oratio = (1.0d * readsetSizes[j]) / ((double) sumOfReadSetSizeOfSupportedPath);
				    if(oratio > 0.7){
					if(HLA.DEBUG3)
					    HLA.log.appendln("oRatio:" + oratio + "\t>\t0.7");
					CustomHashMap oMap = this.paths.get(j).getReadSet();
					IntIterator itr = tMap.keySet().iterator();
					int rmCount=0;
					while(itr.hasNext()){
					    int ci = itr.nextInt();
					    if(oMap.containsKey(ci) || oMap.containsKey(0-ci)){
						if(HLA.DEBUG3)
						    HLA.log.appendln("(" + ci +") in dominant path[" + j+"]" );
						rmCount++;
					    }
					}
					if(rmCount > 0){
					    double updatedRatio = ((readsetSizes[i]-rmCount)*1.0d / ((double) sumOfReadSetSizeOfSupportedPath));
					    if( updatedRatio <= 0.1d ){
						if(HLA.DEBUG3)
						    HLA.log.appendln("Removing due to low updatedRatio:\t" + updatedRatio );
						removalList.add(i);
					    }
					}
				    }
				}
			    }
			}
		    }
		}
	    }
	}
	
	/* Liklihood calculation retain function  --> geared towards retaining the best */
	if(this.paths.size() >= 1){
	    if(HLA.DEBUG3)
		HLA.log.appendln("RUNNING MaxLikelihoodCalcFor BubblePaths");
	    //double[] scoreAndIndices = this.takeMaximumLikeliPair(g);/* [Homo 0-2, Hetero 3-5] 0:homoScore 1:homoIndex1 2:homoIndex2 3:heteroScore 4:heteroIndex1 5: heteroIndex2*/
	    scores = this.takeMaximumLikeliPair(g);
	    double homoScore  = scores.getMaxHomoScore();//scoreAndIndices[0];
	    int homoIndex1 = scores.getMaxHomoGenotypeIndex();//(int) scoreAndIndices[1];
	    int homoIndex2 = homoIndex1;//(int) scoreAndIndices[2];
	    
	    double heteroScore = scores.getMaxHeteroScore();//scoreAndIndices[3];
	    int heteroIndex1 = scores.getMaxHeteroGenotypeIndex1();//(int) scoreAndIndices[4];
	    int heteroIndex2 = scores.getMaxHeteroGenotypeIndex2();//(int) scoreAndIndices[5];
	    int doubleCount1 = scores.getDoubleCountH1();//(int) scoreAndIndices[6];
	    int doubleCount2 = scores.getDoubleCountH2();//(int) scoreAndIndices[7];
	    if(HLA.DEBUG3){
		HLA.log.appendln("2x|H1| = " + doubleCount1);
		HLA.log.appendln("2x|H2| = " + doubleCount2);
	    }
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
		(minCount <= 4 && frequency < 0.15) ||
		(minCount <= 6 && frequency < 0.14) ||
		(minCount <= 8 && frequency < 0.1) ){
		//		(minCount  > 8 && frequency < 0.05) ){
		if(doubleCount1 > doubleCount2)
		    isAlleleWeak2 = true;
		else
		    isAlleleWeak1 = true;
	    }
	    
	    //double[][] genotypeScores; //i: H1 path index, j: H2 path index
	    if(this.paths.size() > 3){
		int[] obHeteroPair = this.getObviousHeteroPair();
		if(obHeteroPair != null){
		    if(obHeteroPair[0] != heteroIndex1 || obHeteroPair[1] != heteroIndex2){
			if(HLA.DEBUG3){
			    HLA.log.appendln("Swapping best heterozygous genotype from [" + heteroIndex1 + ":" + heteroIndex2 
					     + "] --> [" + obHeteroPair[0] + ":" + obHeteroPair[1] + "]" );
			}
			heteroIndex1 = obHeteroPair[0];
			heteroIndex2 = obHeteroPair[1];
			heteroScore = scores.getLogScore(heteroIndex1, heteroIndex2);
			scores.setHeteroIndicies(obHeteroPair[0], obHeteroPair[1]);
			isAlleleWeak1 = false;
			isAlleleWeak2 = false;
		    }
		}
	    }
	    if(HLA.DEBUG3){
		HLA.log.appendln("Best homozygous genotype is :\t" + homoIndex1 + "/" + homoIndex2 + "\tscore:" + homoScore);
		HLA.log.appendln("Best heterozygous genotype is :\t" + heteroIndex1 + "/" + heteroIndex2 + "\tscore:" + heteroScore);
		if(isAlleleWeak1)
		    HLA.log.appendln("H1 count is low : minCount/maxCount" + (minCount/2.0d) + "/" + ((minCount+maxCount)/2.0d));
		else if(isAlleleWeak2)
		    HLA.log.appendln("H2 count is low : minCount/maxCount" + (minCount/2.0d) + "/" + ((minCount+maxCount)/2.0d));
	    }

	    

	    if(HLA.READ_LENGTH >= 200){
		/*
		if(isAlleleWeak1)
		    removalList.add(heteroIndex1);
		else if(isAlleleWeak2)
		    removalList.add(heteroIndex2);
		
		for(int i=0; i<readsetSizes.length;i++){
		    if( (i != heteroIndex1) && (i != heteroIndex2) ){
			removalList.add(i);
			HLA.log.appendln("[Possibly erroneous path] Removing\tPath" + i + "\t");
			this.paths.get(i).printPath();
		    }
		    }*/
		

		double diff = homoScore - heteroScore;
		if(hg.isClassI() && diff > 3){
		    for(int i = 0; i<readsetSizes.length;i++){
			if( i != homoIndex1){
			    removalList.add(i);
			    if(HLA.DEBUG3){
				HLA.log.appendln("[Homo>200][Possibly erroneous path] Removing\tPath" + i + "\t");
				this.paths.get(i).printPath();
			    }
			}
		    }
		}else{
		    if(isAlleleWeak1){
			removalList.add(heteroIndex1);
			if(HLA.DEBUG3){
			    HLA.log.appendln("[Hetero>200Weak1][Possibly erroneous path] Removing\tPath" + heteroIndex1 + "\t");
			    this.paths.get(heteroIndex1).printPath();
			}
		    }else if(isAlleleWeak2){
			removalList.add(heteroIndex2);
			if(HLA.DEBUG3){
			    HLA.log.appendln("[Hetero>200Weak2][Possibly erroneous path] Removing\tPath" + heteroIndex2 + "\t");
			    this.paths.get(heteroIndex2).printPath();
			}
		    }
		    for(int i=0; i<readsetSizes.length;i++){
			if( (i != heteroIndex1) && (i != heteroIndex2) ){
			    removalList.add(i);
			    if(HLA.DEBUG3){
				HLA.log.appendln("[Hetero>200][Possibly erroneous path] Removing\tPath" + i + "\t");
				this.paths.get(i).printPath();
			    }
			}
		    }
		}
	    }else{
		double diff = homoScore - heteroScore;
		if(heteroIndex1 != heteroIndex2){
		    double hetero1Fraction = (readsetSizes[heteroIndex1]) / ((double) sumOfReadSetSizeOfSupportedPath);
		    double hetero2Fraction = (readsetSizes[heteroIndex2]) / ((double) sumOfReadSetSizeOfSupportedPath);
		    int hetIndex1 = heteroIndex1;
		    int hetIndex2 = heteroIndex2;
		    if(hetero1Fraction < hetero2Fraction){
			double tmp = hetero1Fraction;
			hetero1Fraction = hetero2Fraction;
			hetero2Fraction = tmp;
			int tmpI = hetIndex1;
			hetIndex1 = hetIndex2;
			hetIndex2 = tmpI;
		    }
		    if(HLA.DEBUG3){
			HLA.log.appendln("hetero1Fraction:" + hetero1Fraction +"\thetero2Fraction:" + hetero2Fraction);
			HLA.log.appendln("hetIndex1:" + hetIndex1 + "\thetIndex2:" + hetIndex2);
		    }
		    if( (hetero1Fraction + hetero2Fraction) >= 0.8
			&& hetero2Fraction > 0.3){
			for(int i=0; i<readsetSizes.length; i++){
			    if(readsetSizes[i] > 0 && i != hetIndex1 && i != hetIndex2){
				if( (readsetSizes[i] / ((double) sumOfReadSetSizeOfSupportedPath)) < 0.2){
				    removalList.add(i);
				    if(HLA.DEBUG3){
					HLA.log.append("[Possibly errorneous path] Removing\tPath" + i);
					this.paths.get(i).printPath();
				    }
				}
			    }
			}
		    }
		}
		if(HLA.DEBUG3)
		    HLA.log.appendln("diff=homoScore - heteroScore : " + diff);
		if(diff > 0 && diff > 5.3d){ // homoScore > heteroScore
		    for(int i=0; i<readsetSizes.length;i++){
			if(i != homoIndex1){
			    removalList.add(i);
			    if(HLA.DEBUG3){
				HLA.log.appendln("[Possibly erroneous path] Removing\tPath" + i + "\t");
				this.paths.get(i).printPath();
			    }
			}
		    }
		    
		}else{
		    for(int i=0;i<removalList.size();i++){
			if(removalList.getInt(i) == heteroIndex1 || removalList.getInt(i) == heteroIndex2){
			    removalList.removeInt(i);
			    if(HLA.DEBUG3)
				HLA.log.appendln("[Retaing path due to liklihood calculation] Rescueing\tPath" + i + "\t");
			}
		    }
		}
	    }	    

	}else{
	    if(HLA.DEBUG3)
		HLA.log.appendln("NOT Running likelihood calc because there are less than 2 paths remaining.");
	}


	Collections.sort(removalList);

	int pre = -1;
	for(int i=0; i<removalList.size(); i++){
	    if(removalList.getInt(i) == pre){
		removalList.removeInt(i);
		i--;
	    }else
		pre = removalList.getInt(i);
	    
	}

	//removalList is sorted.
	//so we remove from the end so we dont need to worry about index change
	for(int i=removalList.size() - 1; i >= 0; i--){
	    this.paths.get(removalList.getInt(i)).excludePath();
	    this.paths.remove(removalList.getInt(i));
	}
	if(HLA.DEBUG)
	    HLA.log.appendln("Removed (" + removalList.size() + ") paths and\t(" + this.paths.size() + ") left.");
	try{
	    scores.applyRemoval(removalList);
	}catch(Exception e){
	    e.printStackTrace();
	    HLA.log.outToFile();
	    System.exit(-9);
	}
	this.bubbleScores.add(scores);
	return removalList.size();
    }

   
    public PathBaseErrorProb[]  getDataMatrixForLikelihoodCalculation(SimpleDirectedWeightedGraph<Node, CustomWeightedEdge> g){
	PathBaseErrorProb[] pathWisePathBaseErrorProbMatrices = new PathBaseErrorProb[this.paths.size()];
	for(int i=0; i<this.paths.size(); i++)
	    pathWisePathBaseErrorProbMatrices[i] = this.paths.get(i).getBaseErrorProbMatrix(g);
	
	return pathWisePathBaseErrorProbMatrices;
    }
    
    public int[] getObviousHeteroPair(){
	if(this.paths.size() > 2){
	    int max_i = 0;
	    int maxSize = this.paths.get(0).getReadSetSize();
	    int secondMax_i = 1;
	    int secondMaxSize = this.paths.get(1).getReadSetSize();
	    int total = maxSize + secondMaxSize;
	    if(secondMaxSize > maxSize){
		max_i = 1;
		secondMax_i = 0;
		int tmp = secondMaxSize;
		secondMaxSize = maxSize;
		maxSize = tmp;
	    }
	    for(int i=2;i<this.paths.size();i++){
		int curSize = this.paths.get(i).getReadSetSize();
		total += curSize;
		if(curSize > maxSize){
		    secondMax_i = max_i;
		    secondMaxSize = maxSize;
		    max_i = i;
		    maxSize = curSize;
		}else if(curSize > this.paths.get(secondMax_i).getReadSetSize()){
		    secondMax_i = i;
		    secondMaxSize = curSize;
		}
	    }
	    if(total > 30){
		double bestFrac = ( (maxSize + secondMaxSize)*1.0d ) / (total*1.0d);
		double secondFrac = (secondMaxSize*1.0d) / (total*1.0d);
		if(bestFrac > 0.75 && secondFrac > 0.3){
		    int[] result = new int[2];
		    if(max_i < secondMax_i){
			result[0] = max_i;
			result[1] = secondMax_i;
		    }else{
			result[0] = secondMax_i;
			result[1] = max_i;
		    }
		    return result;
		}
	    }
	}
	return null;
    }

    
    //Given all possible paths
    //obtain the pair that gives maxium liklihood of observing data.
    //ONLY considers paths that have phasing-read
    public BubblePathLikelihoodScores takeMaximumLikeliPair(SimpleDirectedWeightedGraph<Node, CustomWeightedEdge> g){
	
	double readFractionScore = 0.0d;
	int readSum = 0;
	
	double curScore = 0.0d;
	//obtain path-wise PathBaseErrorProbMaxtrix
	PathBaseErrorProb[] pathWiseErrorProbMatrices = this.getDataMatrixForLikelihoodCalculation(g);

	for(int i=0; i< pathWiseErrorProbMatrices.length; i++){
	    if(HLA.DEBUG3)
		HLA.log.appendln("path[" + i + "]:\t");
	    char[] pathBases = pathWiseErrorProbMatrices[i].getBases();
	    readSum += pathWiseErrorProbMatrices[i].numReads();
	    if(HLA.DEBUG3){
		for(char c : pathBases){
		    HLA.log.append(c);
		}
		HLA.log.appendln();
	    }
	}

	BubblePathLikelihoodScores scores = new BubblePathLikelihoodScores(this.paths.size());
	
	//getting all possible pairs, including self
	for(int i=0;i<this.paths.size();i++){
	    for(int j=i;j<this.paths.size();j++){
		//NEED TO CALL getScoreForSingleRead over All Reads and sumup the logScore.
		char[] pathBases1 = pathWiseErrorProbMatrices[i].getBases();
		char[] pathBases2 = pathWiseErrorProbMatrices[j].getBases();
		
		//iterate over all reads
		readFractionScore = 0.0d;
		curScore = 0.0d;
		Val whichH = new Val();
		int doubleCountH1 = 0;
		int doubleCountH2 = 0;
		
		double tFraction = (pathWiseErrorProbMatrices[i].numReads() * 1.0d) / (readSum * 1.0d);
		double oFraction = (pathWiseErrorProbMatrices[j].numReads() * 1.0d) / (readSum * 1.0d);
		if(i==j){//homozygous
		    tFraction = tFraction / 2.0d;
		    oFraction = oFraction / 2.0d;
		    readFractionScore += Math.log(tFraction*oFraction);
		}else//heterozygous
		    readFractionScore += Math.log(tFraction*oFraction/HLA.X_FACTOR);
		
		
		for(int k=0; k<this.paths.size();k++){
		    //if( (i== 1 && j == 8) || (i==0 && j == 1) )
		    PathBaseErrorProb curPathErrorMatrix = pathWiseErrorProbMatrices[k];
		    for(int l=0; l<curPathErrorMatrix.numReads(); l++){
			
			boolean debug = false;
			
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
			//if( (i== 0 && j == 1) || (i==0 && j == 2) ){
			
			//HLA.log.appendln("\tP(Di|G) = " + readScore);
			    //			}
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
		
		//HLA.log.appendln("logP( D | Haplotype[ " + i + ":" + j + " ] ) =\t" + curScore);
		if(HLA.DEBUG3)
		    HLA.log.appendln("logP( D | Haplotype[ " + i + ":" + j + " ] ) =\t" + curScore + "\t|H1|x2=" + doubleCountH1 + "\t|H2|x2=" + doubleCountH2);
		scores.updateMax(i, j, curScore, readFractionScore, doubleCountH1, doubleCountH2);
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
	    HLA.log.append("Testing:\t" );
	    for(char c : readBases)
		HLA.log.append(c);
	    HLA.log.append("\nAgainst H1:\t");
	    for(char c : pathBases1)
		HLA.log.append(c);
	    HLA.log.append("\tH2:\t");
	    for(char c : pathBases2)
		HLA.log.append(c);
	    HLA.log.appendln();
	}
	double logProb1 = 0.0d;
	double logProb2 = 0.0d;
	double gapgapmatch = 0.99d;
	double gappenalty = 0.01d;
	
	//for each position in read
	try{
	    for(int i=0; i<readBases.length; i++){
		double errorProb = errorProbs[i];//Math.Pow(0.1, readPhredScores[i]/10.0d);
		double matchProb = 1.0d - errorProb;
		//HLA.log.appendln("read(" + i + "):\terroProb:" + errorProb + "\tmatchProb:" + matchProb );
		double mismatchProb = errorProb / 3.0d;
		char readBase = readBases[i];
		char pathBase1 = pathBases1[i];
		char pathBase2 = pathBases2[i];
		if(debug){
		    HLA.log.appendln("RB:" + readBase  + "\tPB1:" + pathBases1[i] + "\tPB2:" + pathBases2[i]);
		    HLA.log.appendln("MatchP:\t" + matchProb  +"\tMismatchP:\t" + mismatchProb);
		}
		if(readBase == 'N' || pathBase1 == 'N')
		    logProb1 += Math.log(0.25d);//matchProb/4.0d);
		else if(readBase == '-' || pathBase1 == '-'){
		    if(readBase == pathBase1)
			logProb1 += Math.log(gapgapmatch); 
		    else
			logProb1 += Math.log(gappenalty);
		}else if(readBase == pathBase1){
		    if(debug)
			HLA.log.appendln("MatchPB1");
		    logProb1 += Math.log(matchProb);
		}else{
		    if(debug)
			HLA.log.appendln("MismatchPB1");
		    logProb1 += Math.log(mismatchProb);
		}
		if(readBase == 'N' || pathBase2 == 'N')
		    logProb2 += Math.log(0.25d);//matchProb/4.0d);
		else if(readBase == '-' || pathBase2 == '-'){
		    if(readBase == pathBase2)
			logProb2 += Math.log(gapgapmatch); 
		    else
			logProb2 += Math.log(gappenalty);
		}else if(readBase == pathBase2){
		    if(debug)
			HLA.log.appendln("MatchPB2");
		    logProb2 += Math.log(matchProb);
		}else{
		    if(debug)
			HLA.log.appendln("MismatchPB2");
		    logProb2 += Math.log(mismatchProb);
		}
		//HLA.log.appendln("logProb1: " + logProb1 + "\tlogProb2: " + logProb2);
		//logProb += Math.log(avgProb);
		//logProb += columnTransition --> need to add transitionProb
	    }
	}catch(ArrayIndexOutOfBoundsException e){
	    HLA.log.appendln("|eProbs| :" + errorProbs.length);
	    HLA.log.appendln("|rBases| :" + readBases.length);
	    HLA.log.appendln("|pBase1| :" + pathBases1.length);
	    HLA.log.appendln("|pBase2| :" + pathBases2.length);
	    HLA.log.outToFile();
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
	HLA.log.append("logProb1: " + logProb1 + "\tlogProb2: " + logProb2 + "\tsum: " + (logProb1 + logProb2));
	*/
	//HLA.log.appendln("\tlogProb1: " + logProb1 + "\tlogProb2: " + logProb2);
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
	double gapgapmatch = 0.99d;
	double gappenalty = 0.01d;
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
		    logProb1 += Math.log(0.25d);//matchProb/4.0d);
		else if(readBase == '-' || pathBase1 == '-'){
		    if(readBase == pathBase1)
			logProb1 += Math.log(gapgapmatch); 
		    else
			logProb1 += Math.log(gappenalty);
		}else if(readBase == pathBase1)
		    logProb1 += Math.log(matchProb);
		else
		    logProb1 += Math.log(mismatchProb);
		
		if(readBase == 'N' || pathBase2 == 'N')
		    logProb2 += Math.log(0.25d);//matchProb/4.0d);
		else if(readBase == '-' || pathBase2 == '-'){
		    if(readBase == pathBase2)
			logProb2 += Math.log(gapgapmatch); 
		    else
			logProb2 += Math.log(gappenalty);
		}else if(readBase == pathBase2)
		    logProb2 += Math.log(matchProb);
		else
		    logProb2 += Math.log(mismatchProb);
		
		//logProb += Math.log(avgProb);
		//logProb += columnTransition --> need to add transitionProb
	    }
	}catch(ArrayIndexOutOfBoundsException e){
	    HLA.log.appendln("|eProbs| :" + errorProbs.length);
	    HLA.log.appendln("|rBases| :" + readBases.length);
	    HLA.log.appendln("|pBase1| :" + pathBases1.length);
	    HLA.log.appendln("|pBase2| :" + pathBases2.length);
	    e.printStackTrace();
	    System.exit(-1);
	}
	HLA.log.appendln("\tlogProb1: " + logProb1 + "\tlogProb2: " + logProb2);
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
		HLA.log.append("|" + e.getNumActivePath() + "|");
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
	    HLA.log.append("[" + numUnique +  "]");
	}
	for(int i=removalList.size() - 1; i >= 0; i--){
	    this.paths.get(removalList.getInt(i)).excludePath();
	    this.paths.remove(removalList.getInt(i));
	    
	}
	return removalList.size();
    }


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
    
    
    public MergeStatus mergeBubble(Bubble other, int lastSegregationColumnIndex, boolean isClassII, Bubble lastMergedBubble, SimpleDirectedWeightedGraph<Node, CustomWeightedEdge> g){
	
	int pathLength = this.getEnd().get(this.getEnd().size()-1) - this.getStart().get(0);
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
	if(HLA.DEBUG)
	    HLA.log.appendln(">>>>>>>>>>>>> getting intersecti <<<<<<<<<<<<");
	int[][] interBubbleIntersectionSizes = lastMergedBubble.getIntersectionCount(other);

	int[][] interBubbleIntersectionCumulativeSizes = new int[this.paths.size()][other.getPaths().size()]; 

	/* path used counters */
	int[] tpUsed = new int[this.paths.size()];
	int[] opUsed = new int[other.getPaths().size()];

	if(HLA.DEBUG){
	    /* print this paths (DEBUGGIN) */
	    for(int i=0;i<this.paths.size();i++){
		Path tp = this.paths.get(i);
		HLA.log.append("TP(" + i + ")\t<readNum:" + tp.getReadSetSize() + ">\t");
		tp.printInfo();
		tp.printReadSet();
	    }
	    
	    /* print other paths (DEBUGGIN) */
	    for(int i=0;i<other.getPaths().size();i++){
		Path op = other.getPaths().get(i);
		HLA.log.append("OP(" + i + ")\t<readNum:" + op.getReadSetSize() + ">\t");
		op.printInfo();
		op.printReadSet();
	    }
	}
	
	ArrayList<int[]> phasedList = new ArrayList<int[]>();
	ArrayList<Integer> intersectionSizes = new ArrayList<Integer>();
	int intersectionSizesSum = 0;
	int[] intersectionSizesTPSum = new int[this.paths.size()];
	int[] intersectionSizesOPSum = new int[other.getPaths().size()];
	/* check possible paths TP X OP */
	
	for(int i=0;i<this.paths.size();i++){
	    Path tp = this.paths.get(i);
	    for(int j=0; j<other.getPaths().size(); j++){
		Path op = other.getPaths().get(j);
		if(HLA.DEBUG)
		    HLA.log.append("TP(" + i + ")\tX\tOP("+j+"):\t");
		int intersectionSize = tp.isPhasedWith(op);
		interBubbleIntersectionCumulativeSizes[i][j] = intersectionSize;
		if(intersectionSize >= Path.MIN_SUPPORT_PHASING){
		    //if(tp.isPhasedWith(op)){
		    //paths_new.add(tp.mergePaths(op));
		    intersectionSizesSum += intersectionSize;
		    intersectionSizesTPSum[i] += intersectionSize;
		    intersectionSizesOPSum[j] += intersectionSize;
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
		HLA.log.appendln("LOSING TP(" + i + "):\tcurLen:"+pathLength + "\td2ls:" + distanceToLastSegregation);
		if(pathLength >= HLA.READ_LENGTH || distanceToLastSegregation >= (0.5 * HLA.READ_LENGTH)){
		    ms.setSplit(true);
		    HLA.log.appendln("CANT PHASE FURTHER. SPLITTING...");
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
		int lastColumnIndex = this.paths.get(i).getLastUniqueEdgeColumnNumber(g, false);
		int curColumnIndex = other.getStart().get(0);
		int pathSpecificLastSegregation = curColumnIndex - lastColumnIndex + 1;
		//int pathLength = this.getEnd().get(this.getEnd().size()-1) - this.getStart().get(0);
		if(HLA.DEBUG)
		    HLA.log.appendln("LOSING TP(" + i + "):\tcurLen:"+pathLength + "\td2ls:" + distanceToLastSegregation + "\td2psls:"+ pathSpecificLastSegregation);
		if(!otherSignificantSignal){
		    if( (pathLength >= HLA.READ_LENGTH && distanceToLastSegregation >= (0.5 * HLA.READ_LENGTH)) ){
			ms.setSplit(true);
			if(HLA.DEBUG)
			    HLA.log.appendln("[FIRST CHECK]CANT PHASE FURTHER. SPLITTING...");
			//return ms;
		    }else if(pathLength >= HLA.READ_LENGTH && pathSpecificLastSegregation >= (0.75 * HLA.READ_LENGTH)){
			ms.setSplit(true);
			if(HLA.DEBUG)
			    HLA.log.appendln("[LONG SEGREGATION] CANT PHASE FURTHER. SPLITTING...");
		    }else if( distanceToLastSegregation >= (1.5 * HLA.READ_LENGTH) ){
			ms.setSplit(true);
			if(HLA.DEBUG)
			    HLA.log.appendln("[STRONG SIGNAL]CANT PHASE FURTHER. SPLITTING...");
		    }else if(isClassII && pathLength >=200){
			ms.setSplit(true);
			if(HLA.DEBUG)
			    HLA.log.appendln("[CLASS II LENGTH]CANT PHASE FURTHER. SPLITTING...");
		    }else if( this.paths.size() == 2 && phasedList.size() == 1 && pathLength >= 0.5 * HLA.READ_LENGTH && distanceToLastSegregation >= (0.7 * HLA.READ_LENGTH) ){
			ms.setSplit(true);
			if(HLA.DEBUG)
			    HLA.log.appendln("[NECESSARY] CANT PHASE FURTHER. SPLITTING...");
		    }
		}else{
		    if( distanceToLastSegregation >= (1.5 * HLA.READ_LENGTH) ){
			ms.setSplit(true);
			if(HLA.DEBUG)
			    HLA.log.appendln("[STRONG SIGNAL]CANT PHASE FURTHER. SPLITTING...");
		    }else if(isClassII && pathLength >=200){
			ms.setSplit(true);
			if(HLA.DEBUG)
			    HLA.log.appendln("[CLASS II LENGTH]CANT PHASE FURTHER. SPLITTING...");
		    }else if(pathLength >= HLA.READ_LENGTH && pathSpecificLastSegregation >= (0.75 * HLA.READ_LENGTH)){
			ms.setSplit(true);
			if(HLA.DEBUG)
			    HLA.log.appendln("[LONG SEGREGATION] CANT PHASE FURTHER. SPLITTING...");
		    }
		}
	    }
	}
	
	if(ms.isSplit()){
	    if(HLA.DEBUG){
		for(int i=0; i< this.paths.size(); i++){
		    HLA.log.append("TP(" + i + "):\t");
		    this.paths.get(i).printReadSet();
		}
	    }
				 
	    //HLA.log.appendln("CANT PHASE FURTHER. SPLITTING...");
	    return ms;
	}
	/* END OF check for no-phasing paths*/
	if(HLA.DEBUG)
	    HLA.log.appendln("TOTAL of " + phasedList.size() + "\tphased paths.");
	
	int origSizeSum = intersectionSizesSum;

	/*
	for(int j=0; j<opUsed.length; j++){
	    if(opUsed[j] > 1 ){//&& opUsed[j] < phasedList.size()){
		for(int k=0; k<phasedList.size();k++){
		    int[] ijs = phasedList.get(k);
		    if(ijs[1] == j){
			int curSize = intersectionSizes.get(k);
			double d = (intersectionSizes.get(k) * 1.0d) / (origSizeSum * 1.0d);
			HLA.log.appendln("Checking branch:\td:" + d + "\tcurSize:" + curSize + "\tTP(" + ijs[0] + ")\tx\tOP(" + ijs[1] + ")");
			//if( (curSize <2 && d < 0.1 ) || (curSize < 3 && d<0.06) || (curSize >= 3 && d <0.05) ){
			if( (curSize <3 && d < 0.1 ) || (curSize >= 3 && d <0.05) ){  
			    HLA.log.appendln("Pruning branch:\td:"+d+"\tTP(" + ijs[0] + ")\tx\tOP(" + ijs[1] + ")");
			    phasedList.remove(k);
			    intersectionSizesSum -= intersectionSizes.get(k);
			    intersectionSizes.remove(k);
			    k--;//update the index since we took k out.
			    opUsed[j]--;
			    //added 10/04/16
			    tpUsed[ijs[0]]--;
			    if(tpUsed[ijs[0]] == 0){
				int pathLength = this.getEnd().get(this.getEnd().size()-1) - this.getStart().get(0);
				HLA.log.appendln("Pruning results in LOSING TP(" + ijs[0] + "):\tcurLen:"+pathLength + "\td2ls:" + distanceToLastSegregation);
				if(pathLength >= HLA.READ_LENGTH && distanceToLastSegregation >= (0.5 * HLA.READ_LENGTH)){
				    ms.setSplit(true);
				    HLA.log.appendln("CANT PHASE FURTHER. SPLITTING... (PRUNING-INDUCED)");
				    return ms;
				}else if(isClassII && pathLength >=200){
				    ms.setSplit(true);
				    HLA.log.appendln("CANT PHASE FURTHER. SPLITTING... (CLASS II LENGTH)");
				    return ms;
				}
			    }
			    if(tpUsed[ijs[0]] == 1){
				HLA.log.appendln("PRUNING RESULTS IN MORE AGRESSIVE INTERSECTION FOR TP( " + ijs[0] + " )");
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
		double opWiseRatio = (curSize*1.0d) / (intersectionSizesOPSum[ijs[1]] * 1.0d);
		if(HLA.DEBUG)
		    HLA.log.appendln("Checking branch:\td:" + d + "\ttpWiseRatio:" + tpWiseRatio 
				     +  "\topWiseRatio:" + opWiseRatio +  "\tcurSize:" + curSize + "\tTP(" + ijs[0] + ")\tx\tOP(" + ijs[1] + ")"); 
		
		if( ((curSize <3 && d < 0.1 ) || (curSize >= 3 && d <0.05) || (tpWiseRatio < 0.22) 
		     || (curSize >=2 && d <0.15 && origPhasedPathNum > 2 && tpWiseRatio > 0.9 && opWiseRatio < 0.15)
		     || (curSize >=2 && d <0.15 && this.paths.size() == 2 && origPhasedPathNum > 3 && tpWiseRatio < 0.25 && opWiseRatio >=0.65)) ){
		    if( (tpWiseRatio > 0.8 && opWiseRatio > 0.2 && origPhasedPathNum == 2 && pathLength > (HLA.READ_LENGTH/2))
			|| (tpWiseRatio == 1.0d && d >0.075 && pathLength > (HLA.READ_LENGTH*0.7)) ){
			;//dont prune.
		    }else{
			if(HLA.DEBUG)
			    HLA.log.appendln("Pruning branch:\td:"+d+"\tTP(" + ijs[0] + ")\tx\tOP(" + ijs[1] + ")");
			phasedList.remove(i);
			intersectionSizesSum -= curSize;
			intersectionSizes.remove(i);
			i--;
			tpUsed[ijs[0]]--;
			opUsed[ijs[1]]--;
			if(tpUsed[ijs[0]] == 0){
			    int lastColumnIndex = this.paths.get(ijs[0]).getLastUniqueEdgeColumnNumber(g, false);
			    int curColumnIndex = other.getStart().get(0);
			    int pathSpecificLastSegregation = curColumnIndex - lastColumnIndex + 1;
			    if(HLA.DEBUG)
				HLA.log.appendln("lastUniqueColumn:" + lastColumnIndex  + "\tcurColumnIndex:" + curColumnIndex );
			    //int pathLength = this.getEnd().get(this.getEnd().size()-1) - this.getStart().get(0);
			    if(HLA.DEBUG)
				HLA.log.appendln("Pruning results in LOSING TP(" + ijs[0] + "):\tcurLen:"+pathLength + "\td2ls:" + distanceToLastSegregation + "\td2psls:"+ pathSpecificLastSegregation);
			    if(pathLength >= HLA.READ_LENGTH && distanceToLastSegregation >= (0.5 * HLA.READ_LENGTH)){
				ms.setSplit(true);
				if(HLA.DEBUG)
				    HLA.log.appendln("CANT PHASE FURTHER. SPLITTING... (PRUNING-INDUCED)");
				return ms;
			    }else if(pathLength >= HLA.READ_LENGTH && pathSpecificLastSegregation >= (0.6 * HLA.READ_LENGTH)){
				ms.setSplit(true);
				if(HLA.DEBUG)
				    HLA.log.appendln("CANT PHASE FURTHER. SPLITTING... (PRUNING-INDUCED: LONG SEGREGATION)");
				return ms;
			    }else if(isClassII && pathLength >=200){
				ms.setSplit(true);
				if(HLA.DEBUG)
				    HLA.log.appendln("CANT PHASE FURTHER. SPLITTING... (CLASS II LENGTH)");
				return ms;
			    }
			}
			if(HLA.DEBUG){
			    if(tpUsed[ijs[0]] == 1)
				HLA.log.appendln("PRUNING RESULTS IN MORE AGRESSIVE INTERSECTION FOR TP( " + ijs[0] + " )");
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
	    if(tp.getNumUniqueEdges() > 0){
		if(HLA.DEBUG)
		    HLA.log.append("SETTING LAST KNOWN UNIQUE EDGE COLUMN NUMBER AS -->\t");
		tp.setLastKnownUniqueEdgeColumnNumber(tp.getLastUniqueEdgeColumnNumber(g, false));
		if(HLA.DEBUG)
		    HLA.log.appendln(tp.getLastKnownUniqueEdgeColumnNumber());
	    }else{
		if(HLA.DEBUG)
		    HLA.log.appendln("CANT SET BUT LAST KNOWN UNIQUE COLUMN NUMBER IS -->\t" + tp.getLastKnownUniqueEdgeColumnNumber());
	    }
	    //if tp and op are used once, merged path between tp and op is the only PATH
	    if(tpUsed[ijs[0]] == 1 && opUsed[ijs[1]] == 1){
		Path tmpp = tp.mergePathsUnique(op);
		paths_new.add(tmpp);
	    }
	    //tp is used once  op is used multiple times, we pass tp readset.
	    else if(tpUsed[ijs[0]] == 1 && opUsed[ijs[1]] > 1){
		Path tmpp = tp.mergePath1toMany(op);
		paths_new.add(tmpp);
	    }
	    //tp is used multiple time, op is used once, we pass op readset.
	    else if(tpUsed[ijs[0]] > 1 && opUsed[ijs[1]] == 1){
		Path tmpp = tp.mergePathManyto1(op);
		paths_new.add(tmpp);
	    }
	    //tp is used multiple time, op is used multiple times, we pass intersection.
	    else if(tpUsed[ijs[0]] > 1 && opUsed[ijs[1]] > 1){
		Path tmpp = tp.mergePathManytoMany(op);
		paths_new.add(tmpp);
	    }else{
		HLA.log.appendln("SOMETHING IS WRONG. [Bubble.java mergeBubble()]");
		HLA.log.outToFile();
		System.exit(-1);
	    }
	    paths_new.get(paths_new.size() - 1).updateIntersectionSum(intersectionSize, intersectionSizesSum, ijs, interBubbleIntersectionSizes, interBubbleIntersectionCumulativeSizes);
	    
	}
	if(HLA.DEBUG)
	    HLA.log.appendln(paths_new.size() + "\tphased paths in paths_new");
	
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
	    if(HLA.DEBUG){
		HLA.log.appendln("NEW PATHS");
		for(Path p : this.paths)
		    HLA.log.appendln("KnownEdgeCI:\t" + p.getLastKnownUniqueEdgeColumnNumber());
	    }
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
		if(HLA.DEBUG)
		    HLA.log.append("STP(" + i + ")\tX\tSOP("+j+"):\t");
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
