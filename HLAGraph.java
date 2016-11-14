import htsjdk.samtools.SAMRecord;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Set;
import java.util.Arrays;
import java.util.HashSet;
import java.util.Queue;
import java.util.LinkedList;
import java.util.Collection;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;

import org.jgrapht.*;
import org.jgrapht.graph.*;

public class HLAGraph{

    //A B C ... //HLA-DPA1, HLA-DPB1, HLA-DQA1, HLA-DQB1, HLA-DRA, and HLA-DRB1
    private String HLAGeneName; 
    private ArrayList<Sequence> alleles; //
    private HashMap<String, Sequence> alleleHash;
    
    //private SimpleDirectedWeightedGraph<Node, DefaultWeightedEdge> g;
    private SimpleDirectedWeightedGraph<Node, CustomWeightedEdge> g;

    private ArrayList<StringBuffer> interBubbleSequences;
    //private ArrayList<Path> interBubblePaths;

    private ArrayList<TmpPath> interBubblePaths2;

    //private ArrayList<HashMap<Character, Node>> nodeHashList;//list index = columnIndex-1.
    private ArrayList<HashMap<Integer, Node>> nodeHashList;// list index = columnIndex - 1;

    private ArrayList<HLASequence> typingSequences;

    private Node sNode;
    private Node tNode;
    
    private int columnLen;

    private String outputfilename;

    private StringBuffer resultBuffer;

    //keeps track excess lengths added to head and tail of typing regions(exon) due to bubbles in the beginning and end
    private int[] headerExcessLengthBeyondTypingBoundary;
    private int[] tailExcessLengthBeyondTypingBoundary;

    
    /* Outer list index = columnIndex -1 --> insertion point */
    /* Inner list index insertion length */ 
    //private ArrayList<ArrayList<HashMap<Character, Node>>> insertionNodeHashList;
    private ArrayList<ArrayList<HashMap<Integer, Node>>> insertionNodeHashList;

    public void setTypingSequences(ArrayList<HLASequence> seqs){
	this.typingSequences = seqs;
	
	this.writeTypingSequences();
    }

    private void writeTypingSequences(){
	BufferedWriter bw = null;
	try{
	    bw = new BufferedWriter(new FileWriter(HLAGeneName+"_typingSequences_G_group.fa"));
	    for(HLASequence hs : this.typingSequences)
		bw.write(hs.toString());
	    bw.close();
	}catch(IOException ioe){
	    ioe.printStackTrace();
	}
    }

    public SimpleDirectedWeightedGraph<Node, CustomWeightedEdge> getGraph(){
	return this.g;
    }

    /* DO NOT use this to add sNode and tNode */
    public void addVertex(Node n){
	int ci = n.getColIndex();
	this.g.addVertex(n);
	this.nodeHashList.get(n.getColIndex()-1).put(new Integer(n.getIBase()), n);
    }
    
    /* DO NOT use this to add sNode and tNode */
    public void removeVertex(Node n){
	//this.nodeHashList.get(n.getColIndex()-1).remove(new Integer(n.getIBase()));
	this.removeVertexFromNodeHashList(n);
	this.g.removeVertex(n);
    }

    //removes node from nodeHashList. We dont touch insertionNodeHashList
    //because any node added on insertionNodeHashList must have weights.
    private void removeVertexFromNodeHashList(Node n){
	this.nodeHashList.get(n.getColIndex()-1).remove(new Integer(n.getIBase()));
    }

    public void setHLAGeneName(String gn){
	this.HLAGeneName = gn;
    }
    
    
    public Sequence getRefAllele(){
	return this.alleles.get(0);
    }

    public HLAGraph(ArrayList<Sequence> seqs){
	//int numTypingExons = 1;
	//if(this.isClassI())
	//  numTypingExons = 2;
	this.headerExcessLengthBeyondTypingBoundary = new int[2];
	this.tailExcessLengthBeyondTypingBoundary = new int[2];//numTypingExons];
	this.alleles = seqs; 
	this.alleleHash = new HashMap<String, Sequence>();
	for(int i=0;i<this.alleles.size();i++){
	    this.alleleHash.put(this.alleles.get(i).getAlleleName(), this.alleles.get(i));
	}
	//this.g = new SimpleDirectedWeightedGraph<Node, DefaultWeightedEdge>(DefaultWeightedEdge.class);
	this.g = new SimpleDirectedWeightedGraph<Node, CustomWeightedEdge>(CustomWeightedEdge.class);
	this.sNode = new Node('s', 0);
	this.tNode = new Node('t', this.alleles.get(0).getColLength() + 1);
	this.g.addVertex(sNode);
	this.g.addVertex(tNode);
	
	//this.nodeHashList = new ArrayList<HashMap<Character, Node>>();
	this.nodeHashList = new ArrayList<HashMap<Integer, Node>>();
	//this.insertionNodeHashList = new ArrayList<ArrayList<HashMap<Character, Node>>>();
	this.insertionNodeHashList = new ArrayList<ArrayList<HashMap<Integer, Node>>>();
	this.buildGraph();
	this.traverse();
    }

    /*
     * finds all s-t paths in this graph based on BFS technique.
     * Should only be used for each bubble.
     *
     */
    public ArrayList<Path> findAllSTPath(Node s, Node t){
	ArrayList<Path> results = new ArrayList<Path>();
	Queue<Path> pathsQ = new LinkedList<Path>();
	//Set<CustomeWeightedEdge> edges = this.g.outgoingEdgesOf(s);
	Iterator<CustomWeightedEdge> itr = this.g.outgoingEdgesOf(s).iterator();
	//first load all outing edges as paths in paths queue.
	while(itr.hasNext()){
	    pathsQ.add(new Path(itr.next()));
	}
	Path firstPath = null;
	//while we have paths to explore further in the queue
	while((firstPath = pathsQ.poll())!=null){
	    //obtain the vertex at the end for this path
	    Node lastVertex = firstPath.getLastVertex(this.g);
	    //if the last vertex is t, then we add this path in the result
	    if(lastVertex.equals(t)){
		results.add(firstPath);
	    }else{//otherwise, we need to explor the paths further
		itr = this.g.outgoingEdgesOf(lastVertex).iterator();
		while(itr.hasNext()){
		    Path tmpP = firstPath.deepCopy();
		    tmpP.appendEdge(itr.next());
		    pathsQ.add(tmpP);
		}
	    }
	}
	return results;
    }

    public ArrayList<Path> findAllSTPathPruning(Node s, Node t){
	ArrayList<Path> results = new ArrayList<Path>();
	Queue<Path> pathsQ = new LinkedList<Path>();
	Queue<CustomHashMap> readsetQ = new LinkedList<CustomHashMap>(); //we need to keep track of readset to prune branches based on the size
	Iterator<CustomWeightedEdge> itr = this.g.outgoingEdgesOf(s).iterator();
	//first load all outing edges as paths in paths queue.
	while(itr.hasNext()){
	    //Path curP = new Path(itr.next());
	    CustomWeightedEdge curE = itr.next();
	    pathsQ.add(new Path(curE));
	    readsetQ.add(curE.getReadHashSet().clone());
	}
	Path firstPath = null;
	CustomHashMap firstReadSet = null;
	//while we have paths to explore further in the queue
	while((firstPath = pathsQ.poll())!=null){
	    firstReadSet = readsetQ.poll();
	    
	    //obtain the vertex at the end for this path
	    Node lastVertex = firstPath.getLastVertex(this.g);
	    //if the last vertex is t, then we add this path in the result
	    if(lastVertex.equals(t)){
		results.add(firstPath);
	    }else{//otherwise, we need to explor the paths further
		itr = this.g.outgoingEdgesOf(lastVertex).iterator();
		while(itr.hasNext()){
		    Path tmpP = firstPath.deepCopy();
		    CustomHashMap tmpReadSet = firstReadSet.clone();
		    CustomWeightedEdge nextE = itr.next();
		    tmpReadSet.intersectionPE(nextE.getReadHashSet());
		    if(firstReadSet.size() > 0){ // we only add if intersection size is > 0. This greatly prunes paths that are needed to be explored.
			tmpP.appendEdge(nextE);//itr.next());
			pathsQ.add(tmpP);
			readsetQ.add(tmpReadSet);
		    }
		}
	    }
	}
	return results;
    }


    
    //modified so that if pre node is null, create curnode but dont' attempt to connect w/ an edge
    private Node addMissingNode(char b, int colPos, Node cur, Node pre, boolean isRefStrand, byte qual, int readNum){
	cur = new Node(b, colPos);
	
	this.g.addVertex(cur);
	//this.nodeHashList.get(colPos - 1).put(new Character(b), cur);
	this.nodeHashList.get(colPos - 1).put(new Integer(Base.char2ibase(b)), cur);
	if(pre != null){
	    //DefaultWeightedEdge e = this.g.addEdge(pre, cur);
	    this.addAndIncrement(pre, cur, isRefStrand, qual, readNum);
	}//moved readHash to edges
	/*else{
	    //cur.addRead(readNum); 
	    //this.addReadToEdge()
	    }*/
	return cur;
    }
    
    private void addAndIncrement(Node source, Node target, boolean isRefStrand, byte qual, int readNum){
	//target.addRead(readNum); //moved readHash to edges 
	CustomWeightedEdge e = this.g.addEdge(source,target);
	e.addRead(readNum, qual);
	this.g.setEdgeWeight(e, 0.0d);
	e.incrementWeight(this.g, isRefStrand, qual);
    }


    //private void incrementWeight(Node source, Node target){
    private void incrementWeight(Node source, Node target, boolean isRefStrand, byte qual, int readNum){
	//DefaultWeightedEdge e = g.getEdge(source, target);
	//target.addRead(readNum);
	CustomWeightedEdge e = g.getEdge(source, target);
	if(e == null)
	    this.addAndIncrement(source,target, isRefStrand, qual, readNum);
	else{
	    e.addRead(readNum, qual);
	    //target.addRead(readNum); //moved readHash to edges
	    e.incrementWeight(this.g, isRefStrand, qual);//g.setEdgeWeight(e, g.getEdgeWeight(e)+1);
	}
    }
    
    //readNum is a readIdentifier [int]
    public int addWeight(SAMRecord sr, int readNum){
	int numOp = 0;
	Cigar cigar = sr.getCigar();
	byte[] bases = sr.getReadBases(); //ASCII bytes ACGTN=.
	byte[] quals = sr.getBaseQualities();
	for(int i=0; i<quals.length; i++){
	    if(quals[i] < 2)
		quals[i] = 2;
	}
	int baseIndex = 0;
	int refBasePos = sr.getAlignmentStart();
	Node prevnode = null;
	Node curnode = null;
	Base curbase = null;
	Sequence curAllele = this.alleleHash.get(sr.getReferenceName());
	int colPos = curAllele.getColPosFromBasePos(refBasePos);
	boolean isRefStrand = !sr.getReadNegativeStrandFlag();

	/*
	System.err.println(sr.toString());
	System.err.println("start position:\t" + refBasePos);
	System.err.println("Mapped Allele:\t" + sr.getReferenceName());
	System.err.println("Allele Name:\t" + curAllele.getAlleleName());
	System.err.println("CIGAR:\t" + sr.getCigar());
	System.err.println("READ:\t" + sr.getReadString());
	System.err.println("READL:\t" + bases.length);
	System.err.println("ColPos:\t" + colPos);
	for(int i=0; i<bases.length; i++){
	    System.err.print(Base.char2ibase((char)bases[i]));
	}
	System.err.println();
	*/
	//curAllele.printPositions(colPos-1, bases.length);
	
	
	if(cigar==null) return 0;
	for(final CigarElement ce : cigar.getCigarElements()){
	    //System.err.println(ce.toString() + "\t" + ce.getLength());
	    CigarOperator op = ce.getOperator();
	    int cigarLen = ce.getLength();
	    switch(op)
		{
		case S :
		    {
			baseIndex += cigarLen;
			break;
		    }
		case H :
		    break;
		case M :
		    {
			for(int i=0; i<cigarLen; i++){
			    numOp++;
			    /* takes care of jumping over padding area (gaps) in MSA */
			    int tmpColPos = curAllele.getNextColPosForBase(colPos - 1) + 1;
			    if(tmpColPos > colPos){
				for(int j=colPos;j<tmpColPos;j++){
				    HLA.HOPPING++;
				    curnode = this.nodeHashList.get(j-1).get(new Integer(Base.char2ibase('.')));
				    this.incrementWeight(prevnode,curnode,isRefStrand, quals[baseIndex-1], readNum);
				    prevnode=curnode;
				}
				colPos = tmpColPos;
			    }

			    curnode = this.nodeHashList.get(colPos -1).get(new Integer(Base.char2ibase((char)bases[baseIndex])));
			    
			    /* if NO such node is found, we add new node and add edge from prevnode.
			       mismatch that is not covered by reference sequence */
			    if(curnode == null){
				HLA.NEW_NODE_ADDED++;
				curnode = this.addMissingNode((char)bases[baseIndex], colPos, curnode, prevnode, isRefStrand, quals[baseIndex], readNum);
				if(curnode == null)
				    System.err.println("IMPOSSIBLE: curnode NULL again after adding missing node!");
			    }
			    else if(prevnode != null)/* if prevnode is not set. firstBase*/
				this.incrementWeight(prevnode, curnode, isRefStrand, quals[baseIndex], readNum);
			    
			    prevnode=curnode;
			    baseIndex++;
			    //refBasePos++;
			    colPos++;
			    //colPos = curAllele.getNextColPosForBase(colPos - 1) + 1;
			    //colPos++;
			}
			break;
		    }
		case D :
		    {
			for(int i=0; i<cigarLen; i++){
			    numOp++;
			    /* takes care of jumping over padding area (gaps) in MSA */
			    int tmpColPos = curAllele.getNextColPosForBase(colPos - 1) + 1;
			    if(tmpColPos > colPos){
				for(int j=colPos;j<tmpColPos;j++){
				    HLA.HOPPING++;
				    curnode = this.nodeHashList.get(j-1).get(new Integer(Base.char2ibase('.')));
				    this.incrementWeight(prevnode, curnode, isRefStrand, quals[baseIndex-1], readNum);
				    prevnode=curnode;
				}
				colPos = tmpColPos;
			    }
			    /* need to grab gap node at current column */
			    curnode = this.nodeHashList.get(colPos - 1).get(new Integer(Base.char2ibase('.')));
			    
			    /* if NO such node is found, we add new node and add edge from prevnode */
			    if(curnode == null){
				HLA.NEW_NODE_ADDED++;
				curnode = this.addMissingNode('.', colPos, curnode, prevnode, isRefStrand, quals[baseIndex-1], readNum);
			    }else
				this.incrementWeight(prevnode, curnode, isRefStrand, quals[baseIndex-1], readNum);
			    
			    prevnode=curnode;
			    //refBasePos++;
			    colPos++;
			}
			break;
		    }
		case I :
		    {
			// Need to check the colPos distance to nextBase in the allele to see if there are spaces for insertion.
			// If there are spaces for insertion, must insert into those spaces first then insert into insertionNodeHash.
			//
			int insertionIndex = -1;
			for(int i=0; i<cigarLen; i++){
			    HLA.INSERTION++;
			    numOp++;
			    int tmpColPos = curAllele.getNextColPosForBase(colPos - 1) + 1;
			    if(tmpColPos == colPos){//then we must insert into insertionNodeHashList
				insertionIndex++;
				if(this.insertionNodeHashList.get(colPos - 1).size() > insertionIndex){
				    //curnode = this.insertionNodeHashList.get(colPos - 1).get(i).get(new Character((char)bases[baseIndex]));
				    curnode = this.insertionNodeHashList.get(colPos - 1).get(insertionIndex).get(new Integer(Base.char2ibase((char)bases[baseIndex])));
				}else{//we need to add extra position (insertion length)
				    //this.insertionNodeHashList.get(colPos - 1).add(new HashMap<Character, Node>());
				    
				    this.insertionNodeHashList.get(colPos - 1).add(new HashMap<Integer, Node>());
				    curnode = null;
				}
				if(curnode == null){
				    curnode = new Node((char)bases[baseIndex], colPos);
				    HLA.INSERTION_NODE_ADDED++;
				    this.g.addVertex(curnode);
				    this.insertionNodeHashList.get(colPos - 1).get(insertionIndex).put(new Integer(Base.char2ibase((char)bases[baseIndex])), curnode);
				    this.addAndIncrement(prevnode, curnode, isRefStrand, quals[baseIndex], readNum);
				    //DefaultWeightedEdge e = this.g.addEdge(prevnode, curnode);
				    //this.g.setEdgeWeight(e, 0.0d);
				    //this.incrementWeight(prevnode, curnode, isRefStrand,quals[baseIndex]);
				}else{
				    //this.incrementWeight(prevnode, curnode);
				    this.incrementWeight(prevnode, curnode, isRefStrand, quals[baseIndex], readNum);
				}
				prevnode = curnode;
				baseIndex++;
			    }else if(tmpColPos > colPos){//then we must insert here.
				curnode = this.nodeHashList.get(colPos - 1).get(new Integer(Base.char2ibase((char)bases[baseIndex])));
				if(curnode == null){
				    HLA.NEW_NODE_ADDED++;
				    //curnode = this.addMissingNode((char)bases[baseIndex], colPos, curnode, prevnode);
				    curnode = this.addMissingNode((char)bases[baseIndex], colPos, curnode, prevnode, isRefStrand, quals[baseIndex], readNum);
				    if(curnode == null){
					System.err.println("IMPOSSIBLE: curnode NULL again after adding missing node! (1)[addWeight]");
					System.exit(9);
				    }
				}else if(prevnode !=null){
				    HLA.INSERTION_WITH_NO_NEW_NODE++;
				    //this.incrementWeight(prevnode, curnode);
				    this.incrementWeight(prevnode, curnode, isRefStrand, quals[baseIndex], readNum);
				}else if(prevnode == null){
				    System.err.println("SHOULD NOT HAPPEND (2)[addWeight]");//can't start with insertion
				    System.exit(9);
				}

				prevnode = curnode;
				baseIndex++;
				colPos++;
				insertionIndex = -1;
			    }else{//should not happen.
				System.err.println("SHOULD NOT HAPPEND (3)[addWeight]");
				System.exit(9);
			    }

			}
			break;
		    }
		default: System.err.println("UNKNOWN CIGAROP:\t" + ce.toString());
		    break;
		}
	    
	}
	return numOp;
    }
    
    
    private void buildGraph(){
	int numAlleles = this.alleles.size();
	Sequence firstAllele = this.alleles.get(0);
	
	/* for each alleles*/
	//Node sNode = new Node('s', 0);
	//Node tNode = new Node('t', this.alleles.get(0).getColLength() + 1);
	//this.g.addVertex(sNode);
	//this.g.addVertex(tNode);
	for(int i=0; i<numAlleles; i++){
	    //System.err.println("allele " + i);
	    Sequence curSeq = this.alleles.get(i);
	    
	    /* for each base in allele */
	    Node prevNode = sNode;
	    for(int j=0; j<curSeq.getColLength(); j++){
		//System.err.print("[" + j + "]");
		if(i==0){
		    //this.nodeHashList.add(new HashMap<Character, Node>());
		    this.nodeHashList.add(new HashMap<Integer, Node>());
		    
		    //this.insertionNodeHashList.add(new ArrayList<HashMap<Character, Node>>());
		    this.insertionNodeHashList.add(new ArrayList<HashMap<Integer, Node>>());
		}
		//HashMap<Character, Node> curHash = nodeHashList.get(j);
		HashMap<Integer, Node> curHash = nodeHashList.get(j);
		
		//Character curChar = new Character(curSeq.baseAt(j).getBase());
		Integer curInt = new Integer(curSeq.baseAt(j).getIBase());
		
		//Node tmpNode = curHash.get(curChar); //retrieve node
		Node tmpNode = curHash.get(curInt); //retrieve node
		if(tmpNode == null){		//if we have not added this node
		    tmpNode= new Node(curSeq.baseAt(j));
		    this.g.addVertex(tmpNode);
		    //curHash.put(curChar,tmpNode);
		    curHash.put(curInt,tmpNode);
		}
		
		//add an edge
		//DefaultWeightedEdge e;
		CustomWeightedEdge e;
		if(!this.g.containsEdge(prevNode, tmpNode)){
		    e = this.g.addEdge(prevNode,tmpNode);
		    if(prevNode.equals(sNode)){
			//System.err.println("Edge from sNode");
			//if(this.g.getEdge(prevNode,tmpNode) == null)
			//    System.err.println("\tIMPOSSIBLE!!!!");
			//System.err.println("prevNode\t:" + prevNode.toString() + "\tcurNode\t:" + tmpNode.toString());
			this.g.setEdgeWeight(e, Double.MAX_VALUE);
		    }else
			this.g.setEdgeWeight(e, 0.0d);
		}
		prevNode = tmpNode;
	    }
	    //add edge 
	    if(!this.g.containsEdge(prevNode, tNode))
		this.g.setEdgeWeight(this.g.addEdge(prevNode, tNode), Double.MAX_VALUE);
	}
    }

    public void printStartEndNodeInfo(){
	System.err.println(this.sNode.toString() + "|ind("+this.g.inDegreeOf(this.sNode) + ":outd(" + this.g.outDegreeOf(this.sNode ) + ")");
	System.err.println(this.tNode.toString() + "|ind("+this.g.inDegreeOf(this.tNode) + ":outd(" + this.g.outDegreeOf(this.tNode ) + ")");
    }


    public double getTotalWeightForColumn(HashMap<Integer, Node> m, Node preNode){
	double totalWeight = 0;
	Node curNode = null;
	for(int i=0;i<6;i++){
	    curNode = m.get(new Integer(i));
	    CustomWeightedEdge e = this.g.getEdge(preNode,curNode);
	    if(e!=null)
		totalWeight += this.g.getEdgeWeight(e);
	}
	return totalWeight;
    }
    
    public ArrayList<int[]> obtainTypingIntervals(){
	Sequence ref = this.alleles.get(0);
	ArrayList<int[]> typingIntervals = new ArrayList<int[]>();
	if(this.isClassI()){
	    /* typing exon 2 + intron + exon 3 */
	    /*
	    int[] tmp = new int[2];
	    tmp[0] = ref.getBoundaries()[3];
	    tmp[1] = ref.getBoundaries()[6];
	    
	    typingIntervals.add(tmp);
	    */
	    /* typing only exon 2 and 3 */
	    
	    int[] tmp = new int[2];
	    tmp[0] = ref.getBoundaries()[3];
	    tmp[1] = ref.getBoundaries()[4];
	    
	    typingIntervals.add(tmp);
	    
	    int[] tmp2 = new int[2];
	    tmp2[0] = ref.getBoundaries()[5];
	    tmp2[1] = ref.getBoundaries()[6];
	    
	    typingIntervals.add(tmp2);
	    
	}else if (this.isClassII()){
	    int[] tmp2 = new int[2];
	    tmp2[0] = ref.getBoundaries()[3];
	    tmp2[1] = ref.getBoundaries()[4];
	    
	    typingIntervals.add(tmp2);
	}
	return typingIntervals;
    }
    
    public boolean traverseAndWeights(){
	System.err.println("=========================");
	System.err.println("=  " + this.HLAGeneName);
	System.err.println("=========================");

	ArrayList<int[]> typingIntervals = this.obtainTypingIntervals();
	
	Node preNode = null;
	Node curNode = null;
	for(int i=0; i<this.alleles.size(); i++){
	    preNode = null;
	    curNode = null;
	    Sequence curseq = this.alleles.get(i);
	    double exonSum = 0.0d;
	    double exonSump = 0.0d;
	    int exonNumZero = 0;
	    int noEdge = 0;
	    double exonFlow = Double.MAX_VALUE;
	    
	    StringBuffer out = new StringBuffer();
	    out.append(curseq.getAlleleName() + "\n");
	    boolean intact = true;
	    eachallele:
	    for(int j=0; j<typingIntervals.size(); j++){
		int start = typingIntervals.get(j)[0];
		int end = typingIntervals.get(j)[1];
		preNode = null;
		out.append("\nNEXT_EXON\n");
		//need to start a node before the exon start, hence -2, rather than -1 transformation from 1-based to 0-based index
		//k should be 0-based. start and end are 1-based (inclusive, exclusive) index.
		for(int k=start-2; k<end-1; k++){
		    char uchar = Character.toUpperCase(curseq.baseAt(k).getBase());
		    HashMap<Integer, Node> curHash = this.nodeHashList.get(k);
		    curNode = this.nodeHashList.get(k).get(new Integer(Base.char2ibase(uchar)));
		    if(curNode == null){
			preNode = curNode;
			intact = false;
			break eachallele;
		    }
		    if(preNode != null){
			CustomWeightedEdge e = this.g.getEdge(preNode, curNode);
			if(e == null){
			    noEdge++;
			    out.append(uchar + "[NO_EDGE]->");
			    exonFlow = -1.0d;
			    //break;
			}else{
			    double tmpw = this.g.getEdgeWeight(e);
			    double total = this.getTotalWeightForColumn(this.nodeHashList.get(j), preNode);
			    if(tmpw > 0.0d){
				exonSum+=tmpw;
				if(tmpw/total < 0.25d){
				    out.append(("(E)LOWPROB ->\t" + e.getGroupErrorProb() + "\t" + (tmpw/total)) + "\n");
				}else{
				    exonSump+=tmpw;
				}
			    }
			    
			    if(tmpw == 0.0d)
				exonNumZero++;
			    if(tmpw < exonFlow){
				exonFlow = tmpw;
			    }
			    out.append(uchar + "[" + tmpw + "]->");
			}
		    }
		    preNode = curNode;
		}
	    }
	    if(intact){
		out.append(("\n" + curseq.getAlleleName() + "\tNO_EDGE:\t" + noEdge  +"\tE_SUM:\t" + exonSum + "\tE_ZERO:\t" + exonNumZero + "\tE_SUM_P\t" + exonSump + "\tMAXFLOW\t" + exonFlow + "\n"));
	    //out.append(("\n" + curseq.getAlleleName() + "\tSUM:\t" + sum + "\t#ZERO:\t" + numZero + "\tE_SUM:\t" + exonSum + "\tE_ZERO:\t" + exonNumZero + "\tSUM_P:\t" + sump + "\tE_SUM_P\t" + exonSump + "\tMAXFLOW\t" + exonFlow + "\n"));
		System.err.println(out.toString());
	    }
	}
	return true;
    }
    /*
    public boolean traverseAndWeights(){
	System.err.println("=========================");
	System.err.println("=  " + this.HLAGeneName);
	System.err.println("=========================");
	
	Node preNode;
	Node curNode;
	//double exonFlow = Double.MAX_VALUE;
	for(int i=0; i<this.alleles.size(); i++){
	    preNode = this.sNode;
	    Sequence curseq = this.alleles.get(i);
	    double sum = 0.0d;
	    double sump = 0.0d;
	    int numZero = 0;
	    
	    double exonSum = 0.0d;
	    double exonSump = 0.0d;
	    int exonNumZero = 0;
	    double exonFlow = Double.MAX_VALUE;
	    
	    //System.err.println(curseq.getAlleleName());
	    StringBuffer out = new StringBuffer();
	    out.append(curseq.getAlleleName() + "\n");
	    boolean intact = true;
	    for(int j=0; j<curseq.getColLength(); j++){
		char uchar = Character.toUpperCase(curseq.baseAt(j).getBase());
		HashMap<Integer, Node> curHash = this.nodeHashList.get(j);
		curNode = this.nodeHashList.get(j).get(new Integer(Base.char2ibase(uchar)));
		if(!preNode.equals(this.sNode)){
		    //System.err.print(uchar + "[" + this.g.getEdgeWeight(this.g.getEdge(preNode, curNode)) + "]->");
		    CustomWeightedEdge e = this.g.getEdge(preNode, curNode);
		    if(e == null){
			intact = false;
			break;
		    }
		    double tmpw = this.g.getEdgeWeight(this.g.getEdge(preNode, curNode));
		    double total = this.getTotalWeightForColumn(this.nodeHashList.get(j), preNode);
		    if(tmpw > 0){
			if(tmpw/total < 0.3d){
			    ;//System.err.println("(I)LOWPROB ->\t" + this.g.getEdge(preNode,curNode).getGroupErrorProb()+ "\t" + (tmpw/total));
			}else{
			    sump+=tmpw;
			}
		    }
		    sum+=tmpw;
		    if(curseq.withinTypingExon(j+1)){
			if(tmpw == 0.0d)
			    exonNumZero++;
			if(tmpw < exonFlow){
			    //System.err.print("*FU*");
			    exonFlow = tmpw;
			}
			exonSum+=tmpw;
			if(tmpw > 0){
			    if(tmpw/total < 0.3d){
				out.append(("(E)LOWPROB ->\t" + this.g.getEdge(preNode,curNode).getGroupErrorProb() + "\t" + (tmpw/total)) + "\n");
				//System.err.println("(E)LOWPROB ->\t" + this.g.getEdge(preNode,curNode).getGroupErrorProb() + "\t" + (tmpw/total));
			    }else{
				exonSump+=tmpw;
			    } 
			}
			out.append(uchar + "[" + this.g.getEdgeWeight(this.g.getEdge(preNode, curNode)) + "]->");//System.err.print(uchar + "[" + this.g.getEdgeWeight(this.g.getEdge(preNode, curNode)) + "]->");
		    }
		    if(tmpw == 0.0d){
			numZero++;
			//if(curseq.withinTypingExon(j+1)){
			//  exonNumZero++;
			//}
		    }
		    
		}
		preNode = curNode;
	    }
	    if(intact){
		out.append(("\n" + curseq.getAlleleName() + "\tSUM:\t" + sum + "\t#ZERO:\t" + numZero + "\tE_SUM:\t" + exonSum + "\tE_ZERO:\t" + exonNumZero + "\tSUM_P:\t" + sump + "\tE_SUM_P\t" + exonSump + "\tMAXFLOW\t" + exonFlow + "\n"));
		System.err.println(out.toString());
	    }
	    //System.err.println("\n" + curseq.getAlleleName() + "\tSUM:\t" + sum + "\t#ZERO:\t" + numZero + "\tE_SUM:\t" + exonSum + "\tE_ZERO:\t" + exonNumZero + "\tSUM_P:\t" + sump + "\tE_SUM_P\t" + exonSump + "\tMAXFLOW\t" + exonFlow);
	}
	return true;
    }
    */    
    
    public void traverse(){
	System.err.println("Traversing (" + this.alleles.size() + ")");
	Node preNode;// = this.sNode;
	Node curNode;
	for(int i=0; i<this.alleles.size(); i++){
	    this.alleles.get(i).verify();
	    preNode = this.sNode;
	    Sequence curseq = this.alleles.get(i);
	    for(int j=0; j<curseq.getColLength(); j++){
		//System.err.println("Traversing [" + i + "," + j + "]");
		char uchar = Character.toUpperCase(curseq.baseAt(j).getBase());
		char lchar = Character.toUpperCase(curseq.baseAt(j).getBase());
		
		//HashMap<Character, Node> curHash = this.nodeHashList.get(j);
		HashMap<Integer, Node> curHash = this.nodeHashList.get(j);
		
		//if(curHash.get(new Character(uchar)) != null){
		if(curHash.get(new Integer(Base.char2ibase(uchar))) != null){
		    //System.err.println("NODE FOUND IN HASH[UPPER}");
		    //curNode = curHash.get(new Character(uchar));
		    curNode = curHash.get(new Integer(Base.char2ibase(uchar)));
		    /*
		    if(this.g.getEdge(preNode, curNode) == null)
			System.err.println("\tWRONG, THIS SHOULD ALREADY BE IN THE GRAPH.\n" + "prevNode\t:" + preNode.toString() + "\tcurNode\t:" + curNode.toString());
		    else
			System.err.println("Weight : " + this.g.getEdgeWeight(this.g.getEdge(preNode,curNode)));
		    */
		    preNode = curNode;
			
		    //}else if(curHash.get(new Character(lchar)) != null){
		}else if(curHash.get(new Integer(Base.char2ibase(lchar))) != null){
		    //System.err.println("NODE FOUND IN LOWER}");
		    //curNode = curHash.get(new Character(lchar));
		    curNode = curHash.get(new Integer(Base.char2ibase(lchar)));
		    /*
		    if(this.g.getEdge(preNode, curNode) == null)
			System.err.println("\tWRONG, THIS SHOULD ALREADY BE IN THE GRAPH.");
		    else
			System.err.println("Weight : " + this.g.getEdgeWeight(this.g.getEdge(preNode,curNode)));
		    */
		    preNode = curNode;
		}else{
		    ;//System.err.println("NODE NOT FOUND IN THH GRAPH");
		}
	    }
	}
	System.err.println("DONE Traversing");
    }
    
    public void updateEdgeWeightProb(){
	Set<CustomWeightedEdge> eSet = g.edgeSet();
	Iterator<CustomWeightedEdge> itr = eSet.iterator();
	CustomWeightedEdge e = null;
	while(itr.hasNext()){
	    e = itr.next();
	    e.computeGroupErrorProb();
	    //System.err.println(e.toString());
	}
    }

    public boolean isClassI(){
	if( this.HLAGeneName.equals("A") 
	    || this.HLAGeneName.equals("B") 
	    || this.HLAGeneName.equals("C") 
	    ){
	    return true;
	}
	return false;
    }
    
    public boolean isClassII(){
	if( //this.HLAGeneName.equals("DPA1") 
	    //|| this.HLAGeneName.equals("DPB1") 
	    this.HLAGeneName.equals("DQA1") 
	    || this.HLAGeneName.equals("DQB1") 
	    //|| this.HLAGeneName.equals("DRA") 
	    || this.HLAGeneName.equals("DRB1") 
	    ){
	    return true;
	}
	return false;
    }
    

    /*
     * countBubbles() --> returns ArrayList of simple bubbles (NOT merged)
     *                    and sets interbubbleSequences in this class.
     *
     * processBubbles() --> merges bubbles by checking reads support information.
     */
    public void countBubblesAndMerge(StringBuffer rb){
	this.resultBuffer = rb;
	this.processBubbles(this.countBubbles());
    }

    /*
    public void processBubblesOLD(ArrayList<Bubble> bubbles){

	for(int i=0; i<bubbles.size(); i++){
	    bubbles.get(i).initBubbleSequences();
	}

	Bubble superBubble = bubbles.get(0);
	
	for(int i=0; i<bubbles.size(); i++){
	    System.err.println("B" + i + "|");
	    bubbles.get(i).printPaths();
	}
	
	superBubble.printBubbleSequence();
	
	System.err.println("(iteration 0):\t" + superBubble.getNumPaths());
	
	for(int i=1; i<bubbles.size(); i++){
	    System.err.println("\t(attempting merging)\t" + bubbles.get(i).getNumPaths());
	    bubbles.get(i).printBubbleSequence();
	    System.err.print("(SB)\t");
	    superBubble.printBubbleSequenceSizes(); 
	    System.err.print("(OB)\t");
	    bubbles.get(i).printBubbleSequenceSizes();
	    superBubble.mergeBubble(bubbles.get(i));
	    System.err.println("**********************************");
	    superBubble.printBubbleSequenceSizes();
	    System.err.println("**********************************");
	    superBubble.printBubbleSequence();
	    System.err.println("(iteration " + i + "):\t" + superBubble.getNumPaths());
	}
	
	superBubble.printResults(this.interBubbleSequences);
    }
    */
    public void processBubbles(ArrayList<Bubble> bubbles){
	/* to load actual bubble sequence in each paths found in each bubble */
	System.err.println("**************************");
	System.err.println("Checking numBubbles: " + bubbles.size());
	for(int i=0; i<bubbles.size(); i++){
	    if(bubbles.get(i).isFirstBubble()){
		System.err.println("Bubble (" + i + "):\t[FB]" );
	    }
	    bubbles.get(i).initBubbleSequences();
	}
	
	/* superBubble is a merged bubbles. Ideally, you want to have just one bubble. */
	ArrayList<Bubble> superBubbles = new ArrayList<Bubble>();
	
	
	Bubble curSuperBubble = bubbles.get(0);
	Bubble lastMergedBubble = curSuperBubble;
	int lastSegregationColumnIndex = curSuperBubble.getStart().get(0);

	System.err.println("(iteration 0):\t" + curSuperBubble.getNumPaths());
	
	for(int i=1; i<bubbles.size(); i++){
	    System.err.println("\t(attempting merging)\t" + bubbles.get(i).getNumPaths());
	    bubbles.get(i).printBubbleSequence();
	    System.err.print("(SB)\t");
	    curSuperBubble.printBubbleSequenceSizes(); 
	    System.err.print("(OB)\t");
	    bubbles.get(i).printBubbleSequenceSizes();
	    //boolean phased = curSuperBubble.mergeBubble(bubbles.get(i));
	    MergeStatus ms = null;
	    if(!bubbles.get(i).isFirstBubble()){
		ms = curSuperBubble.mergeBubble(bubbles.get(i), lastSegregationColumnIndex, this.isClassII(), lastMergedBubble);
		lastMergedBubble = bubbles.get(i);
	    }
	    //if we are cutting here
	    if(bubbles.get(i).isFirstBubble() || ms.isSplit()){
		if(bubbles.get(i).isFirstBubble())
		    System.out.println("NOT PHASING OVER DIFFERENT EXONS --> setting OB as curSuperBubble");
		else
		    System.out.println("CANT PHASE --> setting OB as curSuperBubble.");
		superBubbles.add(curSuperBubble);
		curSuperBubble = bubbles.get(i);
		lastMergedBubble = curSuperBubble;
		//need to update segregationColumnIndex
		lastSegregationColumnIndex = curSuperBubble.getStart().get(0);
	    }
	    //if not cutting
	    else{
		//if we have a segreation, need to updated segregationColumnIndex
		if(ms.isSegregating())
		    lastSegregationColumnIndex = ms.getLastSegregationColumnIndex();
		
		System.err.println("**********************************");
		curSuperBubble.printBubbleSequenceSizes();
		System.err.println("**********************************");
		curSuperBubble.printBubbleSequence();
	    }
	    /*
	    if(!phased){
		System.out.println("NOT PHASED --> setting OB as curSuperBubble.");
		superBubbles.add(curSuperBubble);
		curSuperBubble = bubbles.get(i);
	    }else{
		System.err.println("**********************************");
		curSuperBubble.printBubbleSequenceSizes();
		System.err.println("**********************************");
		curSuperBubble.printBubbleSequence();
		}*/
	    System.err.println("(iteration " + i + "):\t" + curSuperBubble.getNumPaths());
	}
	
	superBubbles.add(curSuperBubble);

	System.err.println("\n\n<---------------------------------->\nCHECKING INTER-SUPERBUBBLE PHASING:\n<---------------------------------->\n");
	this.checkSuperBubbleLinkages(superBubbles);
	
	//this.printBubbleResults(superBubbles, bubbles);
	//this.compareInterBubbles(superBubbles);
	ArrayList<ArrayList<AllelePath>> fracturedPaths = this.getFracturedPaths(superBubbles, bubbles);
	
	this.allelePathPrintTest(fracturedPaths);//print test of fractured candidate. print super bubble sequences
	this.allelePathToFastaFile(fracturedPaths);//writes superbubble sequences as fasta file
	
	ArrayList<SuperAllelePath> superpaths = this.generateSuperAllelePaths(fracturedPaths); 
	this.superAllelePathToFastaFile(superpaths); //writes full length candidate allele concatenating super bubbles as fasta file
	this.printScoreForMaxLikeliPair(superpaths, superBubbles);
	this.pathAlign(superpaths); // aligns to DB for typing.
	
    }


    public void printScoreForMaxLikeliPair(ArrayList<SuperAllelePath> superpaths, ArrayList<Bubble> superBubbles){
	//allProduct, jointProduct, avgProduct, MAXFLOW
	double[] curBest = {Double.NEGATIVE_INFINITY, Double.NEGATIVE_INFINITY, Double.NEGATIVE_INFINITY, 0.0d};
	double[] curSecondBest = {Double.NEGATIVE_INFINITY, Double.NEGATIVE_INFINITY, Double.NEGATIVE_INFINITY, 0.0d};
	
	int[][] bestIndicies = new int[4][2];
	int[][] secondBestIndicies = new int[4][2];
	for(int i = 0; i<superpaths.size(); i++){
	    for(int j=i; j<superpaths.size(); j++){
		double[] scores = superpaths.get(i).getJointProbability(superpaths.get(j), superBubbles);
		double[] jointWeightFlow = superpaths.get(i).jointTraverse(superpaths.get(j), this.g);
		System.err.println("AllelePair [" + i + ":" + j + "]\t{ + " 
				   + scores[0] + "\t" 
				   + scores[1] + "\t" 
				   + scores[2] + "\t" 
				   + scores[3] + "\t" 
				   + scores[4] + "\t" 
				   + scores[5] 
				   + "\tE_SUM:" + jointWeightFlow[0] 
				   + "\tMAXFLOW:" + jointWeightFlow[1]
				   + "}");
		//higher the better 
		for(int k=0; k<3; k++){
		    if(curBest[k] < scores[k+3]){
			curSecondBest[k] = curBest[k];
			curBest[k] = scores[k+3];
			secondBestIndicies[k][0] = bestIndicies[k][0];
			secondBestIndicies[k][1] = bestIndicies[k][1];
			bestIndicies[k][0] = i;
			bestIndicies[k][1] = j;
		    }else if(curSecondBest[k] < scores[k+3]){
			curBest[k] = scores[k+3];
			bestIndicies[k][0] = i;
			bestIndicies[k][1] = j;
		    }
		}
		if(curBest[3] < jointWeightFlow[1]){
		    curSecondBest[3] = curBest[3];
		    curBest[3] = jointWeightFlow[1];
		    secondBestIndicies[3][0] = bestIndicies[3][0];
		    secondBestIndicies[3][1] = bestIndicies[3][1];
		    bestIndicies[3][0] = i;
		    bestIndicies[3][1] = j;
		}else if(curSecondBest[3] < jointWeightFlow[1]){
		    curBest[3] = jointWeightFlow[1];
		    bestIndicies[3][0] = i;
		    bestIndicies[3][1] = j;
		}
		
		
	    }
	}
	System.err.println("-------- AllProductMetric --------");
	System.err.print("RANK 1:\t");
	this.printBest(bestIndicies, curBest, 0);
	System.err.println("RANK 2:\t");
	this.printBest(secondBestIndicies, curSecondBest, 0);
	
	System.err.println("-------- JointProductMetric --------");
	System.err.print("RANK 1:\t");
	this.printBest(bestIndicies, curBest, 1);
	System.err.println("RANK 2:\t");
	this.printBest(secondBestIndicies, curSecondBest, 1);
	
	System.err.println("-------- AvgProductMetric --------");
	System.err.print("RANK 1:\t");
	this.printBest(bestIndicies, curBest, 2);
	System.err.println("RANK 2:\t");
	this.printBest(secondBestIndicies, curSecondBest, 2);

	System.err.println("-------- JointMaxFlowMetric --------");
	System.err.print("RANK 1:\t");
	this.printBest(bestIndicies, curBest, 3);
	System.err.println("RANK 2:\t");
	this.printBest(secondBestIndicies, curSecondBest, 3);
	
	/* superAllelePath-wise best score printing */
	/*
	  int count = 0;
	for(SuperAllelePath sap : superpaths){
	    double[] weightFlow = sap.traverse(this.g);
	    System.err.println("Superpath[" + count + "]\tE_SUM:" + weightFlow[0] + "\tMAXFLOW:" + weightFlow[1]);
	    count++;
	    }*/
    }

    private void printBest(int[][] indicies, double[] curBest, int typeIndex){
	System.err.println("AllelePari[" + indicies[typeIndex][0] + ":" + indicies[typeIndex][1] + "]\t{AP:"
			   + curBest[0] + "\tJP:" + curBest[1] + "\tAVP:" + curBest[2] + "\tMF" + curBest[3] + "}"
			   );
    }
    
    
    public void pathAlign(ArrayList<SuperAllelePath> superpaths){
	int count = 1;
	
	for(SuperAllelePath sap : superpaths){
	    String candidate = sap.getSequenceBuffer().toString();//p.toString(this.g, count);//, this.headerExcessLengthBeyondTypingBoundary, this.tailExcessLengthBeyondTypingBoundary);
	    count++;
	    String sapname = sap.toSimpleString();
	    String subject = null;
	    //String maxName = null;
	    ArrayList<String> maxName = new ArrayList<String>();
	    int maxIdenticalLen = 0;
	    //ArrayList<Integer> maxIdenticalLen = new ArrayList<Integer>();
	    //Result maxR = null;
	    ArrayList<Result> maxR =new ArrayList<Result>();
	    
	    boolean foundPerfect = false;
	    for(HLASequence subjscan : this.typingSequences){
		subject = subjscan.getSequence();
		if(candidate.equals(subject)){
		    Result curR = new Result(candidate.length(), subject);
		    //maxIdenticalLen = curR.getIdenticalLen();
		    //maxName = subj.getGroup().getGroupString();
		    maxR.add(curR);
		    maxName.add(subjscan.getGroup().getGroupString());
		    System.err.println("Found perfect match.");
		    foundPerfect = true;
		    break;
		}
	    }
	    
	    if(!foundPerfect){
		for(HLASequence subj : this.typingSequences){
		    subject = subj.getSequence();
		    Result curR = NWAlign.runDefault(candidate, subject);
		    /*if(subj.getGroup().getGroupString().equals("A*01:01:01G")){
		      System.err.println(candidate);
		      System.err.println(subject);
		      System.err.println("A*01:01:01G\t" + curR.toString());
		      }*/
		    if(curR.getIdenticalLen() >= maxIdenticalLen){
			if(curR.getIdenticalLen() > maxIdenticalLen){
			    maxName = new ArrayList<String>();
			    maxIdenticalLen = curR.getIdenticalLen();
			    maxR =new ArrayList<Result>();
			}
			//maxName.add(subj.getGroup().getGroupString());
			//maxIdenticalLen.add(curR.getIdenticalLen());
			maxName.add(subj.getGroup().getGroupString());
			maxR.add(curR);
		    }
		    
		}
	    }
	    //System.err.print("BEST MATCH:" + );
	    for(int i=0;i<maxR.size();i++){
		System.err.println("["+ sapname+  "]BEST MATCH:\t" + maxName.get(i) + "\t" + maxR.get(i).getIdenticalLen() + "\t" + maxR.get(i).getIdentity());
		this.resultBuffer.append(maxName.get(i) + "\t" + maxR.get(i).getIdenticalLen() + "\t" + maxR.get(i).getIdentity() + "\t" + maxR.get(i).getScore() + sapname + "\n");
	    }
	    //System.err.println("BEST MATCH:\t" + maxName + "\t" + maxIdenticalLen + "\t" + maxR.getIdentity());
	    //this.resultBuffer.append(maxName + "\t" + maxIdenticalLen + "\t" + maxR.getIdentity() + "\t" + maxR.getScore() + sapname+"\n");
	    
	    //this.resultBuffer.append(maxR.toAlignmentString() + "\n");
	}
    }

    /*
    public void printBubbleResults(ArrayList<Bubble> superBubbles){
	int startIndex = 0;

	System.out.println("Printing\t" + superBubbles.size() + "\tfractured super bubbles.");
	
	int count = 0;
	for(Bubble sb : superBubbles){
	    System.out.println("\tSuperBubble\t" + count);
	    startIndex = sb.printResults(this.interBubbleSequences, startIndex);
	    count++;
	}

	
    }
    */


    public void checkSuperBubbleLinkages(ArrayList<Bubble> superBubbles){
	ArrayList<int[]>[] pLists = new ArrayList[(superBubbles.size()-1)*superBubbles.size()/2];
	int count = 0;
	/* for each superBubble*/
	for(int i=0;i<superBubbles.size();i++){
	    Bubble sb_i = superBubbles.get(i);
	    /* pairing with another superBubble */
	    for(int j=i+1; j<superBubbles.size();j++){
		Bubble sb_j = superBubbles.get(j);
		//int[0]: path index for first bubble
		//int[1]: path index for second bubble
		//int[2]: number of reads supporting this phasing path
		ArrayList<int[]> phasedList = sb_i.getPhasedSuperBubbles(sb_j);
		pLists[count] = phasedList;
		count++;
		if(phasedList.size() > 0){
		    System.err.println("Phasing evidence FOUND between SB(" + i + ") : SB(" + j + ")" );
		    for(int[] index : phasedList)
			System.err.println("SB(" + i + ")-" + index[0] + " : SB(" + j + ")-" + index[1]);
		}else
		    System.err.println("NO phasing evidence between SB(" + i + ") : SB(" + j + ")" );
		
	    }
	}
    }

    public void selectGreedyForSuperBubbleLinking(ArrayList<int[]>[] phasedLists){
	//for(ArrayList<int[]>)
    }

    /*
    public void compareInterBubbles(ArrayList<Bubble> superBubbles){
	//System.out.println(">>>>>>>>>>>>>>>> Checking interbubbles  <<<<<<<<<<");
	//for(int i=0; i<this.interBubbleSequences.size();i++){
	//    System.out.println("[I" + i + "]:\t" + this.interBubbleSequences.get(i).toString() + "\t" + this.interBubblePaths.get(i).toSimplePathString(this));
	//    }

	int k = 0;
	for(int i=0; i<superBubbles.size(); i++){
	    Bubble sb = superBubbles.get(i);
	    Path firstPath = sb.getPaths().get(0);
	    ArrayList<StringBuffer> bubbleSequences = firstPath.getBubbleSequences();
	    ArrayList<CustomWeightedEdge> orderedEdgeList = firstPath.getOrderedEdgeList();
	    int curEdgePos = 0;
	    int curMaxPos = 0;
	    for(int j=0; j<bubbleSequences.size(); j++){
		System.out.println("[I" + k + "]:\t" + this.interBubbleSequences.get(k).toString() + "\t" + this.interBubblePaths.get(k).toSimplePathString(this));
		k++;
		System.out.print("[B:" + j +"]" + bubbleSequences.get(j).toString() + "\t");
		curMaxPos += sb.getBubbleLengths().get(j).intValue();
		for(;curEdgePos < curMaxPos; curEdgePos++){
		    System.out.print(this.g.getEdgeTarget(orderedEdgeList.get(curEdgePos)).getBase());
		}
		System.out.println();
	    }
	}
    }
    */
    public void setFileName(String f){
	this.outputfilename = f;
    }
    
    public ArrayList<DNAString> generateCandidates(ArrayList<ArrayList<DNAString>> fracturedSequences){
	
	ArrayList<DNAString> sequences = new ArrayList<DNAString>();
	for(DNAString ds : fracturedSequences.get(0)){
	    sequences.add(ds.deepCopy());
	}
	
	//for superBubble
	for(int i=1; i<fracturedSequences.size(); i++){
	    ArrayList<DNAString> otherSequences = fracturedSequences.get(i);
	    ArrayList<DNAString> results = new ArrayList<DNAString>();
	    for(int j=0; j < sequences.size(); j++){
		for(int k=0; k < otherSequences.size(); k++){
		    results.add(sequences.get(j).mergeDeep(otherSequences.get(k)));
		}
	    }
	    sequences = results;
	}

	BufferedWriter bw = null;
	try{
	    bw = new BufferedWriter(new FileWriter(this.outputfilename + "_" + this.HLAGeneName + ".typed.fa.candidates"));
	    for(DNAString seq : sequences)
		bw.write(seq.toFasta().toString());
	    bw.close();
	}catch(IOException ioe){
	    ioe.printStackTrace();
	}   
	
	return sequences;
    }
    
    
    /*
    public ArrayList<ArrayList<Path>> mergePathsOverSuperBubbles(ArrayList<Bubble> superBubbles){
	int startIndex = 0;
	int count = 0;
	ArrayList<ArrayList<Path>> fracturedPaths = new ArrayList<ArrayList<Path>>();
	
	for(Bubble sb : superBubbles){
	    ArrayList<Path> paths = new ArrayList<Path>();
	    fracturedPaths.add(paths);
	    
	}
    }
    */
    /*
    public void printBubbleResults(ArrayList<Bubble> superBubbles, ArrayList<Bubble> bubbles){
	//StringBuffer output = new StringBuffer();
	int startIndex = 0;
	
	System.out.println("Printing\t" + superBubbles.size() + "\tfractured super bubbles.");
	//output.append(superBubbles.size() + "\tfractured SuperBubbles\n");
	int count = 0;

	
	//over each super bubble
	ArrayList<ArrayList<DNAString>> fracturedSequences = new ArrayList<ArrayList<DNAString>>();

	int bubbleOffset = 0;
	Bubble pre = null;
	for(Bubble sb : superBubbles){
	    if(pre != null){
		bubbleOffset += pre.numBubbles();
	    }
	    ArrayList<DNAString> sequences = new ArrayList<DNAString>();
	    fracturedSequences.add(sequences);
	    System.out.println("\tSuperBubble\t" + count);
	    System.out.println("\t\tbubbleOffset:\t" + bubbleOffset);
	    startIndex = sb.printResults(this.interBubbleSequences, startIndex, sequences, this.HLAGeneName , count,  bubbles, bubbleOffset);
	    count++;
	    pre = sb;
	}

	BufferedWriter bw = null;
	try{
	    bw = new BufferedWriter(new FileWriter(this.outputfilename + "_" + this.HLAGeneName + ".typed.fa"));
	    for(ArrayList<DNAString> fseq : fracturedSequences){
		for(DNAString ds : fseq)
		    bw.write(ds.toFasta().toString());
	    }
	    //bw.write(output.toString());
	    bw.close();
	}catch(IOException ioe){
	    ioe.printStackTrace();
	}

	ArrayList<DNAString> candidateAlleles = this.generateCandidates(fracturedSequences);
	//this.candidateAlign(candidateAlleles);
    }
    */

    public ArrayList<ArrayList<AllelePath>> getFracturedPaths(ArrayList<Bubble> superBubbles, ArrayList<Bubble> bubbles){
	int startIndex = 0;
	int count = 0;

	System.out.println("Printing\t" + superBubbles.size() + "\tfractured super bubbles.");
	//inner list holds paths found for one superBubble
	//outer list holds multiple superBubbles
	ArrayList<ArrayList<AllelePath>> fracturedPaths = new ArrayList<ArrayList<AllelePath>>();

	int bubbleOffset = 0;
	Bubble presb = null;
	int sbIndex = 0;
	for(Bubble sb : superBubbles){
	    if(presb != null){
		bubbleOffset += presb.numBubbles();
	    }
	    ArrayList<AllelePath> paths = new ArrayList<AllelePath>();
	    fracturedPaths.add(paths);
	    //NEED TO ADD TRIM FUNCTIONALITY FOR HEADER AND TAIL BUBBLES!!! --> Trim function ADDED
	    startIndex = sb.mergePathsInSuperBubbles(this.interBubblePaths2, startIndex, paths, this.HLAGeneName, count, this.g, bubbles, bubbleOffset);
	    count++;
	    presb = sb;
	}
	
	return fracturedPaths;
	//this.pathPrintTest(this.generateCandidatePaths(fracturedPaths));
	//this.pathAlign(this.generateCandidatePaths(fracturedPaths));
    }
    
    
    public ArrayList<SuperAllelePath> generateSuperAllelePaths(ArrayList<ArrayList<AllelePath>> fracturedSequences){
	ArrayList<SuperAllelePath> superpaths = new ArrayList<SuperAllelePath>();
	
	//for(AllelePath ap : fracturedSequences.get(0))
	for(int i=0; i<fracturedSequences.get(0).size();i++){
	    AllelePath ap = fracturedSequences.get(0).get(i);
	    superpaths.add(new SuperAllelePath(this.HLAGeneName));
	    superpaths.get(superpaths.size()-1).addAllelePath(ap, i);
	}
	for(int i=1; i<fracturedSequences.size(); i++){
	    ArrayList<AllelePath> nextSequences = fracturedSequences.get(i);
	    ArrayList<SuperAllelePath> results = new ArrayList<SuperAllelePath>();
	    for(int j=0; j<superpaths.size(); j++){
		for(int k=0; k < nextSequences.size(); k++){
		    results.add(superpaths.get(j).clone());
		    results.get(results.size()-1).addAllelePath(nextSequences.get(k), k);
		}
	    }
	    superpaths = results;
	}
	
	return superpaths;
    }

    public void allelePathPrintTest(ArrayList<ArrayList<AllelePath>> fracturedAllelePaths){
	for(int i=0; i<fracturedAllelePaths.size(); i++){
	    ArrayList<AllelePath> paths = fracturedAllelePaths.get(i);
	    System.out.println("SUPER BUBBLE [" + i + "]");
	    for(int j=0; j<paths.size(); j++){
		AllelePath ap = paths.get(j);
		ap.printPath(this.g, i, j);
	    }
	}
    }
    
    public void allelePathToFastaFile(ArrayList<ArrayList<AllelePath>> fracturedAllelePaths){
	BufferedWriter bw = null;
	try{
	    bw = new BufferedWriter(new FileWriter(this.outputfilename + "_" + this.HLAGeneName + ".typed.fa"));
	    for(ArrayList<AllelePath> faps : fracturedAllelePaths){
		for(AllelePath ap : faps){
		    bw.write(ap.toFasta().toString());
		}
		//bw.close();
	    }
	    bw.close();
	}catch(IOException ioe){
	    ioe.printStackTrace();
	}
    }

    public void superAllelePathToFastaFile(ArrayList<SuperAllelePath> superAllelePaths){
	BufferedWriter bw = null;
	try{
	    bw = new BufferedWriter(new FileWriter(this.outputfilename + "_" + this.HLAGeneName + ".typed.fa.candiates"));
	    for(SuperAllelePath sap : superAllelePaths)
		bw.write(sap.toFasta().toString());
	    bw.close();
	}catch(IOException ioe){
	    ioe.printStackTrace();
	}
    }
    

    
    /*
    public void getFracturedPathsOLD(ArrayList<Bubble> superBubbles, int[] headerExcessArr, int[] tailExcessArr){
	int startIndex = 0;
	int count = 0;
	
	//inner list holds paths found for one superBubble
	//outer list holds multiple superBubbles
	ArrayList<ArrayList<Path>> fracturedPaths = new ArrayList<ArrayList<Path>>();
	Bubble presb = null;
	ArrayList<Path> prePaths = null;
	Bubble sb = null;
	int firstBubbleCount = 0;
	int headerExcess,tailExcess;
	//for(sb : superBubbles){
	for(int i=0;i<superBubbles.size(); i++){
	    sb = superBubbles.get(i);
	    ArrayList<Path> paths = new ArrayList<Path>();
	    fracturedPaths.add(paths);
	    startIndex = sb.mergePathsInSuperBubbles(this.interBubblePaths, startIndex, paths, this.HLAGeneName, count);
	    if(sb.isFirstBubble()){
		headerExcess = headerExcessArr[firstBubbleCount];
		tailExcess = (firstBubbleCount > 0 ? tailExcessArr[firstBubbleCount-1] : 0);
		if(presb != null){
		    for(Path p : prePaths)
			p.trimExcess(0, tailExcess);
		}
		for(Path p: paths)
		    p.trimExcess(headerExcess, 0);
		firstBubbleCount++;
		presb = sb;
		prePaths = paths;
	    }

	    count++;
	}
	if(sb !=null && tailExcessArr[firstBubbleCount-1] > 0){
	    for(Path p : prePaths)
		p.trimExcess(0, tailExcessArr[firstBubbleCount-1]);
	}
	
	//this.pathPrintTest(this.generateCandidatePaths(fracturedPaths));
	this.pathAlign(this.generateCandidatePaths(fracturedPaths));
    }
    

    
    public void pathPrintTest(ArrayList<Path> ps){
	int count = 1;
	for(Path p : ps){
	    p.printPath(this.g, count);//, this.headerExcessLengthBeyondTypingBoundary, this.tailExcessLengthBeyondTypingBoundary);
	    count++;
	}
	}
    */
    
    public void candidateAlign(ArrayList<DNAString> candidates){
	int count = 1;
	for(DNAString candidateDNA : candidates){
	    String candidate = candidateDNA.getSequence();
	    String subject = null;
	    String maxName = null;
	    String maxHit = null;
	    int maxIdenticalLen = 0;
	    Result maxR = null;
	    for(HLASequence subj : this.typingSequences){
		subject = subj.getSequence();
		Result curR = NWAlign.runDefault(candidate, subject);
		/*if(subj.getGroup().getGroupString().equals("A*01:01:01G")){
		    System.err.println(candidate);
		    System.err.println(subject);
		    System.err.println("A*01:01:01G\t" + curR.toString());
		    }*/
		if(curR.getIdenticalLen() >= maxIdenticalLen){
		    maxIdenticalLen = curR.getIdenticalLen();
		    maxName = subj.getGroup().getGroupString();
		    maxR = curR;
		    maxHit = subject;//curR.getHit();
		    if(curR.getIdentity() == 1.0d){
			System.err.println("Found perfect match.");
			break;
		    }
		}
	    }

	    System.err.println("BEST MATCH:\t" + maxName + "\t" + maxIdenticalLen + "\t" + maxR.getIdentity());
	    System.err.println("Query:\n"+candidate);
	    System.err.println("Hit:\n"+maxHit);
	    
	    this.resultBuffer.append(maxName + "\t" + maxIdenticalLen + "\t" + maxR.getIdentity() + "\t" + maxR.getScore() + "\n");
	    
	    this.resultBuffer.append(maxR.toAlignmentString() + "\n");
	}
    }

    /*
    public void pathAlign(ArrayList<Path> ps){
	int count = 1;
	for(Path p : ps){
	    String candidate = p.toString(this.g, count);//, this.headerExcessLengthBeyondTypingBoundary, this.tailExcessLengthBeyondTypingBoundary);
	    count++;
	    String subject = null;
	    String maxName = null;
	    int maxIdenticalLen = 0;
	    Result maxR = null;
	    for(HLASequence subj : this.typingSequences){
		subject = subj.getSequence();
		Result curR = NWAlign.runDefault(candidate, subject);
		if(curR.getIdenticalLen() >= maxIdenticalLen){
		    maxIdenticalLen = curR.getIdenticalLen();
		    maxName = subj.getGroup().getGroupString();
		    maxR = curR;
		    if(curR.getIdentity() == 1.0d){
			System.err.println("Found perfect match.");
			break;
		    }
		}
	    }

	    System.err.println("BEST MATCH:\t" + maxName + "\t" + maxIdenticalLen + "\t" + maxR.getIdentity());
	    this.resultBuffer.append(maxName + "\t" + maxIdenticalLen + "\t" + maxR.getIdentity() + "\t" + maxR.getScore() + "\n");
	    
	    this.resultBuffer.append(maxR.toAlignmentString() + "\n");
	}
	
    }


    public ArrayList<Path> generateCandidatePaths(ArrayList<ArrayList<Path>> fracturedPaths){
	ArrayList<Path> paths = new ArrayList<Path>();
	//add paths of the first superBubble
	for(Path p : fracturedPaths.get(0)){
	    paths.add(p.deepCopy());
	}
	
	//for each of next superBubble
	for(int i=1; i<fracturedPaths.size(); i++){
	    ArrayList<Path> otherPaths = fracturedPaths.get(i);
	    ArrayList<Path> results = new ArrayList<Path>();
	    //for each current path
	    for(int j=0; j < paths.size(); j++){
		//for each next option
		for(int k=0; k < otherPaths.size(); k++){
		    results.add(paths.get(j).combinePaths(otherPaths.get(k)));
		}
	    }
	    paths = results;
	}
	return paths;
    }
    
    
    public void selectBestHits(ArrayList<DNAString> candidates){
	ArrayList<Integer> score = new ArrayList<Integer>();
	for(DNAString seq:candidates){
	    score.add(findBestHit(seq));
	}
    }

    public int findBestHit(DNAString seq){
    	int score = 0;
	//run alignment
	return score;
    }
    */
    public ArrayList<Bubble> countBubbles(){
	System.err.println("=========================");
	System.err.println("=  " + this.HLAGeneName);
	System.err.println("=========================");

	ArrayList<Bubble> bubbles = new ArrayList<Bubble>();

	ArrayList<int[]> typingIntervals = this.obtainTypingIntervals();

	/* counters */
	int numBubbles = 0;
	int curBubbleLength = 1;
	int lastStartOfBubble = 0;
	//ArrayList<Integer> numPaths = new ArrayList<Integer>();
	ArrayList<Integer> bubbleLengths = new ArrayList<Integer>(); // keeps track of bubble lengths. Bubble length is length excluding collapsing nodes. L-2
	ArrayList<Integer> coordinates = new ArrayList<Integer>(); //keeps track of start coordinates of bubbles
	/* counters */
	
	Node curSNode = null;

	this.interBubbleSequences = new ArrayList<StringBuffer>();
	//this.interBubblePaths = new ArrayList<Path>();
	
	this.interBubblePaths2  = new ArrayList<TmpPath>();

	StringBuffer curbf = new StringBuffer("");
	TmpPath tp = new TmpPath();
	for(int i=0; i<typingIntervals.size(); i++){
	    int start = typingIntervals.get(i)[0];
	    int end = typingIntervals.get(i)[1];
	    
	    curBubbleLength = 1;
	    lastStartOfBubble = start - 2;
	    //boolean headerBubble = false;

	    boolean firstBubble = true; // to demarcate the first bubble of the interval

	    //Node preNode = null;

	    int k;
	    /* FOR EACH POSITION in a TYPING INTERVAL*/
	    for(k=start-1;k<end-1;k++){
		HashMap<Integer, Node> columnHash = this.nodeHashList.get(k);
		Integer[] keys = columnHash.keySet().toArray(new Integer[0]);
		
		/*it's a collapsing node if curBubbleLength > 2
		  else it's a possible start of bubble.*/
		if(keys.length == 1){
		    //headerBubble = false;
		    /* then it must be a collapsing node; */
		    if(curBubbleLength > 1){
			this.interBubbleSequences.add(curbf);
			//this.interBubblePaths.add(tp.toPath(this.g));
			this.interBubblePaths2.add(tp);
			//this.interBubblePaths.add(curP);
			curBubbleLength++;
			numBubbles++;
			//numPaths.add(new Integer(this.analyzeBubble(lastStartOfBubble, k)));
			bubbleLengths.add(new Integer(curBubbleLength-2));
			coordinates.add(new Integer(lastStartOfBubble));
			if(firstBubble){
			    //if(i>0)//if it's not first interval, we need to update last bubble
				//	bubbles.get(bubbles.size()-1).trimPaths(0,this.tailExcessLengthBeyondTypingBoundary[i-1]);
			    bubbles.add(new Bubble(this, curSNode, columnHash.get(keys[0]), firstBubble, this.headerExcessLengthBeyondTypingBoundary[i], 0));
			    //bubbles.get(bubbles.size()-1).trimPath(this.headerExcessLengthBeyongTypingBoundary[i], 0);
			    firstBubble = false;
			}else
			    bubbles.add(new Bubble(this, curSNode, columnHash.get(keys[0])));
			curSNode = columnHash.get(keys[0]);
			//preNode = curSNode;
			lastStartOfBubble = k;
			curBubbleLength = 1;
			//curP = new Path();
			curbf = new StringBuffer("");
			curbf.append(curSNode.getBase());
			tp = new TmpPath();
			tp.appendNode(curSNode);
		    }
		    /* Possible Start of a Bubble or straight path */
		    else{
			curSNode = columnHash.get(keys[0]);
			curbf.append(curSNode.getBase());
			tp.appendNode(curSNode); 
			/*if(prNode == null)
			    preNode = curSNode;
			else{
			    curP.appendEdge(this.g.getEdge(preNode, curSNode));
			    preNode = curSNode;
			    }*/
			lastStartOfBubble = k;
			curBubbleLength = 1;
		    }
		}else if(keys.length > 1){//middle of bubble

		    /* NEED TO FIX THIS TO ALLOW BUBBLE TO BE USED at the boundaries*/
		    if(k==(start-1)){// || headerBubble){
			System.err.println("[k] = " + k);
			int tmpBubbleLength = 1;
			for(int l=start-2;;l--){
			    System.err.println("trying new k: [k] = " + l);
			    tmpBubbleLength++;
			    HashMap<Integer, Node> tmpHash = this.nodeHashList.get(l);
			    Integer[] tmpKeys = tmpHash.keySet().toArray(new Integer[0]);
			    if(tmpKeys.length == 1){
				System.err.println("Found the new start!");
				curSNode = tmpHash.get(tmpKeys[0]);
				curbf.append(curSNode.getBase());// this is actually unecessary
				//curbf=new StringBuffer("");
				tp.appendNode(curSNode);
				lastStartOfBubble = l;
				curBubbleLength = tmpBubbleLength;
				this.headerExcessLengthBeyondTypingBoundary[i] = curBubbleLength - 1;
				System.err.println("Setting Trimming length(header):\t" + this.headerExcessLengthBeyondTypingBoundary[i]);
				break;
			    }
			}
			//this.interBubbleSequences.add(new StringBuffer(""));
			//headerBubble = true;
			/*curSNode = columnHash.get(keys[0]);
			curbf.append(curSNode.getBase());
			tp.appendNode(curSNode);
			lastStartOfBubble = k;
			curBubbleLength = 1;
			*/
		    }else{ //mid-bubble: just increment bubble length
			curBubbleLength++;
			//preNode = null;
		    }
		}else{//disconnected graph.
		    System.err.println("This should NOT HAPPEN");
		}
	    }
	    //need to update here to handle "End-Bubble" (bubble sitting at the end and not concluded)
	    if(curBubbleLength > 1){
		System.err.println(">>>>>>>Bubble at the end:\t[curBubbleLength]:"+ curBubbleLength);
		int preLength = curBubbleLength;
		for(;;k++){
		    HashMap<Integer, Node> columnHash = this.nodeHashList.get(k);
		    Integer[] keys = columnHash.keySet().toArray(new Integer[0]);
		    curBubbleLength++;
		    if(keys.length == 1){
			this.interBubbleSequences.add(curbf);
			//this.interBubblePaths.add(tp.toPath(this.g));
			this.interBubblePaths2.add(tp);
			System.err.println("Found the new end!");
			numBubbles++;
			bubbleLengths.add(new Integer(curBubbleLength-2));
			coordinates.add(new Integer(lastStartOfBubble));
			//if(firstBubble){
			//   bubbles.add(new Bubble(this, curSNode, columnHash.get(keys[0]), firstBubble));
			//    firstBubble = false;
			//}else
			this.tailExcessLengthBeyondTypingBoundary[i] = curBubbleLength - preLength;
			System.err.println("Setting Trimming length(tail):\t" + this.tailExcessLengthBeyondTypingBoundary[i]);
			bubbles.add(new Bubble(this, curSNode, columnHash.get(keys[0]), false, 0, this.tailExcessLengthBeyondTypingBoundary[i]));
			curSNode = columnHash.get(keys[0]);
			lastStartOfBubble = k;
			curBubbleLength = 1;
			curbf = new StringBuffer("");
			curbf.append(curSNode.getBase());
			tp = new TmpPath();
			tp.appendNode(curSNode);
			
			break;
		    }
		}
	    }//else{
	    this.interBubbleSequences.add(curbf);
	    //this.interBubblePaths.add(tp.toPath(this.g));
	    this.interBubblePaths2.add(tp);
	    
	    curbf = new StringBuffer("");
	    tp = new TmpPath();
	//
	    /*
	    this.interBubbleSequences.add(curbf);
	    this.interBubblePaths.add(tp.toPath(this.g));
	    curbf = new StringBuffer("");
	    tp = new TmpPath();
	    if(curBubbleLength > 1){
		System.err.println(">>>>>>>Bubble at the end:\t[curBubbleLength]:"+ curBubbleLength);
	    }
	    */
	}
	System.err.println("NumBubbles:\t" + numBubbles + "\tfound");
	for(int i=0; i<bubbleLengths.size(); i++){
	    System.err.print(bubbleLengths.get(i).intValue() + "\t");
	}
	System.err.println();
	for(int i=0; i<bubbleLengths.size(); i++){
	    System.err.print(coordinates.get(i).intValue() + "\t");
	}
	System.err.println();
	
	return bubbles;
    }

    

    //write code to find number of paths and 
    //return the number of paths in the bubble.
    //move column-wise and update number of paths going through each vertex.
    /*
    private int analyzeBubble(int start, int end){
	
	Integer[] keys = this.nodeHashList.get(start).keySet().toArray(new Integer[0]);
	
	
	for(int i=start+1; i<=end; i++){
	    //HashMap<Integer, Node> columnHash = this.nodeHashList.get(i);
	    //Integer[] keys = columnHash.keySet().toArray(new Integer[0]);
	    this.updateNumPathFwd(i-1, i);
	}
	return 0;
	}*/

    //update numPathFwd in current column
    /*
    private void updateNumPathFwd(int pre, int cur){
	Collection<Node> preNodes = this.nodeHashList.get(pre).values();
	Collection<Node> curNodes = this.nodeHashList.get(cur).values();
	
	Iterator<Node> curItr = curNodes.iterator();
	while(curItr.hasNext()){
	    Node curNode = curItr.next();
	    Iterator<Node> preItr = preNodes.iterator();
	    while(preItr.hasNext()){
		Node preNode = preItr.next();
		if(this.g.getEdge(preNode, curNode) != null){
		    curNode.incrementNumPathInBubbleFwd(preNode.getNumInBubbleFwd());
		}
	    }
	}
    }
    */

    
    public void countBubbles(boolean typingExonOnly){
	int startIndex, endIndex;
		
	if(typingExonOnly){
	    int[] boundaries = this.alleles.get(0).getBoundaries();
	    if(this.alleles.get(0).isClassI()){//if class I : type exon 2 and 3
		startIndex = boundaries[3];
		endIndex = boundaries[6];
	    }else{// if class II : type exon 2
		startIndex = boundaries[3];
		endIndex = boundaries[4];
	    }
	}else{
	    startIndex = 0;
	    endIndex = this.nodeHashList.size();
	}
	
	int numBubbles = 0;

	Node sNode = new Node(4, startIndex);
	ArrayList<Node> preNodes = new ArrayList<Node>();
	preNodes.add(sNode);
	boolean preStart = true;
	
	
	int bubbleSize = 1;
	int numPath = 1;
	
	for(int i = startIndex; i <= endIndex; i++){
	    
	    HashMap<Integer, Node> curHash = this.nodeHashList.get(i);
	    //Set<Integer> keyset = curHash.keySet();
	    Integer[] keys = curHash.keySet().toArray(new Integer[0]);      
	    
	    if(keys.length == 1){//only one option --> it's collaping node or part of just a straight path
		if(bubbleSize > 1){//if bublleSize > 1, then it's the end end of bubble
		    numBubbles++;     
		    System.err.println("Bubble[" + numBubbles + "]:Size(" + bubbleSize + "):numPath(" + numPath + ")" );
		    preNodes = new ArrayList<Node>();
		    preNodes.add(curHash.get(keys[0]));
		    preStart = false;
		    bubbleSize = 1;
		    numPath = 1;
		}else{
		    preNodes = new ArrayList<Node>();
		    preStart = false;
		}
	    }else if(keys.length > 1){
		//checking previous column nodes to this column node
		for(int p=0; p < preNodes.size(); p++){
		    Node pNode = preNodes.get(p);
		    int branching=0;
		    for(int q=0; q<keys.length; q++){
			Node qNode = curHash.get(keys[q]);
			CustomWeightedEdge e = this.g.getEdge(pNode, qNode);
			if(e != null && this.g.getEdgeWeight(e) > 0)
			    branching++;
		    }
		    if(branching > 2){
			if(preStart){
			    numPath += (branching - 1);
			}else{
			    int ind = this.g.inDegreeOf(pNode);
			    numPath += ind*branching - ind;
			}
		    }
		}
	    }
	}
    }

    

    //insertionNodes are indexed at same position as endColumns
    //meaning: insertionNodes should be inserted in between startColumns and endColumns.
    public void flattenInsertionNodes(){
	
	ArrayList<int[]> typingIntervals = this.obtainTypingIntervals();

	int fCount = 0;	
	
	for(int i=typingIntervals.size()-1; i>-1; i--){
	    int start = typingIntervals.get(i)[0];
	    int end   = typingIntervals.get(i)[1];
	    
	    for(int j=end-1; j >= start; j--){
		int insSize = this.insertionNodeHashList.get(j).size();
		//there is insertion, we need to flatten.
		if(insSize > 0 && this.isThereConnectionToInsertionNodes(insSize, j)){
		    fCount++;
		    this.shiftColumnsByInsertionSize(insSize, j);
		}
	    }
	}
	
	System.err.println(this.HLAGeneName + "\t>>>>> FLATTENED InsertionBubble:\t" + fCount );
    }

    
    //fromColumnIndex is 0-based columnIndex
    private boolean isThereConnectionToInsertionNodes(int insSize, int fromColumnIndex){
	System.err.println("[isThereConnection] Checking at fromColumnIndex : " + fromColumnIndex + "\tInsSize: " + insSize);

	HashMap<Integer, Node> startNodes = nodeHashList.get(fromColumnIndex-1);
	boolean sConnection = false;
	boolean eConnection = false;
	HashMap<Integer, Node> sInsHash = this.insertionNodeHashList.get(fromColumnIndex).get(0);
	HashMap<Integer, Node> eInsHash = this.insertionNodeHashList.get(fromColumnIndex).get(insSize - 1);
	HashMap<Integer, Node> endNodes = nodeHashList.get(fromColumnIndex);
	
	System.out.println("[isThereConnectionToInsertionNodes] HashIndex: " + (fromColumnIndex - 1) );
	sConnection = this.isThereConnection(startNodes, sInsHash);
	eConnection = this.isThereConnection(eInsHash, endNodes);
	if(sConnection || eConnection){
	    if(sConnection)
		System.err.println("[isThereConnection] connection between startNodes and sInsHash found!");
	    else
		System.err.println("[isThereConnection] NO connection between startNodes and sInsHash found!");
	    if(eConnection)
		System.err.println("[isThereConnection] connection between eInsHash and endNodes found!");
	    else
		System.err.println("[isThereConnection] NO connection between eInsHash and endNodes found!");
	}
	return sConnection && eConnection;
    }
    
    //just to check if there edges between s and t
    private boolean isThereConnection(HashMap<Integer, Node> s, HashMap<Integer, Node> t){
	Integer[] sKeys = new Integer[0];
	sKeys = s.keySet().toArray(sKeys);
	Integer[] eKeys = new Integer[0];
	eKeys = t.keySet().toArray(eKeys);
	
	for(int i=0;i<sKeys.length; i++){
	    if(sKeys[i].intValue() != 4){
		for(int j=0; j<eKeys.length; j++){
		    //System.err.print("eKyes[j] intval\t");
		    //System.err.println(eKeys[j].intValue());
		    if(eKeys[j].intValue() != 4){
			//System.out.println("Actual index value in node: " + s.get(sKeys[i]).getColIndex());
			CustomWeightedEdge e = this.g.getEdge(s.get(sKeys[i]), t.get(eKeys[j]));
			if(e != null)
			    return true;
		    }
		}
	    }		
	}
	return false;
    }

    /* fromColumnIndex is 0-based index --> this is where insertion happens */
    /* 0based(List index):          0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 */
    /* 1based(CI in Node and Base): 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 */
    /* from ColumnIndex at 5, insSize of 2*/

    private void shiftColumnsByInsertionSize(int insSize, int fromColumnIndex){
	
	HashMap<Integer, Node> startNodes = nodeHashList.get(fromColumnIndex-1);
	HashMap<Integer, Node> endNodes = nodeHashList.get(fromColumnIndex);
	
	Iterator<Integer> itr_s = startNodes.keySet().iterator();
	System.err.println("\n**STARTNODES:");
	while(itr_s.hasNext()){
	    System.err.println(startNodes.get(itr_s.next()).toString());
	}

	Iterator<Integer> itr_e = endNodes.keySet().iterator();
	System.err.println("\n**ENDNODES:");
	while(itr_e.hasNext()){
	    System.err.println(endNodes.get(itr_e.next()).toString());
	}	

	Node pre = null;
	Node[] gapNodes = new Node[insSize];
	
	/* here we first shift endNodes and all nodes after that by insSize
	 * to acquire insSize many column space for insertionNodeHashList.
	 */
	for(int i=0; i<insSize;i++){
	    HashMap<Integer, Node> insHash_i = this.insertionNodeHashList.get(fromColumnIndex).get(i);
	    this.adjustColumnIndex(insHash_i, fromColumnIndex + i + 1);//1-base column position 
	    nodeHashList.add(fromColumnIndex+i, insHash_i); //insert insHash_i
	    Node cur = new Node('.', fromColumnIndex + i + 1); // 1-base column position
	    this.addVertex(cur);//add vertex and add it to nodeHashList;
	    if(pre !=null)
		this.g.addEdge(pre,cur);
	    gapNodes[i] = cur;
	    pre = cur;
	}
	
	System.err.println("checking edges between gapNodes[]");
	for(int i=1;i<gapNodes.length;i++){
	    CustomWeightedEdge e = this.g.getEdge(gapNodes[i-1], gapNodes[i]);
	    if(e == null)
		System.err.println("No edges found between between gapNodes["+(i-1) + "] and gapNodes[" + i + "]");
	    
	}

	/* adding spaces to Alleles as well */
	for(int i=0; i<this.alleles.size(); i++)
	    this.alleles.get(i).insertBlanks(fromColumnIndex, insSize);
	
	
	/* we shift all columns after insertion, so updating all columnIndex */
	for(int i=fromColumnIndex+insSize; i<this.nodeHashList.size(); i++)
	    this.adjustColumnIndex(this.nodeHashList.get(i), i+1);//need to updated with 1-base column position
	
	/* remove all edges between start node and end nodes and re-route them through gap nodes by adding new edges and assign weights and readset accordingly*/
	double weightSum = this.getWeightSumsBetween2Columns(startNodes, endNodes, gapNodes);

	/* DEBUGGING prints*/
	itr_s = startNodes.keySet().iterator();
	System.err.println("\n**STARTNODES:");
	while(itr_s.hasNext()){
	    System.err.println(startNodes.get(itr_s.next()).toString());
	}

	System.err.println("**CONNECTED NODES TO START-GAP:");
	CustomWeightedEdge[] inEdges = this.g.incomingEdgesOf(gapNodes[0]).toArray(new CustomWeightedEdge[1]);
	for(CustomWeightedEdge e : inEdges){
	    System.err.println(this.g.getEdgeSource(e).toString());
	}

	itr_e = endNodes.keySet().iterator();
	System.err.println("\n**ENDNODES:");
	while(itr_e.hasNext()){
	    System.err.println(endNodes.get(itr_e.next()).toString());
	}

	System.err.println("**CONNECTED NODES TO END-GAP:");
	CustomWeightedEdge[] outEdges = this.g.outgoingEdgesOf(gapNodes[gapNodes.length -1]).toArray(new CustomWeightedEdge[1]);
	for(CustomWeightedEdge e : outEdges){
	    System.err.println(this.g.getEdgeTarget(e).toString());
	}
	
    }

    private void shiftColumnsByInsertionSizeOLD(int insSize, int fromColumnIndex){
	
	HashMap<Integer, Node> startNodes = nodeHashList.get(fromColumnIndex-2);
	HashMap<Integer, Node> endNodes = nodeHashList.get(fromColumnIndex-1);

	//we need to insert <insSize>-many columns first
	Node pre = null;
	Node[] gapNodes = new Node[insSize];
	ArrayList<Base> insBases = new ArrayList<Base>();//new Base[insSize];


	//insert insSize-many columns with gapNodes and transfer insertionNodes to nodeHashList.
	for(int i=0; i<insSize; i++){
	    //add a space first then add the vertex --> gets the space(HashMap) from insertionNodeHashList
	    HashMap<Integer, Node> insHash_i = this.insertionNodeHashList.get(fromColumnIndex-1).get(i);
	    this.adjustColumnIndex(insHash_i, fromColumnIndex + i);//this.adjustColumnIndex(insHash_i, fromColumnIndex + i + 1);
	    nodeHashList.add(fromColumnIndex + i, insHash_i);
	    Node cur = new Node('.', fromColumnIndex + i + 1);
	    this.addVertex(cur);//add vertex and add to nodeHashList
	    if(pre != null)
		this.g.addEdge(pre, cur);
	    gapNodes[i] = cur;
	    pre = cur;
	    insBases.add(new Base('-', 0,0,0,true,1));
	}

	
	/* adding spaces to Alleles as well*/
	for(int i=0; i<this.alleles.size(); i++){
	    //this.alleles.get(i).insertBlanks(fromColumnIndex, insBases);
	    this.alleles.get(i).insertBlanks(fromColumnIndex-1, insSize);
	}
	/*
	//insert insSize-many columns with gapNodes
	for(int i=0; i<insSize; i++){
	    //add a space first then add the vertex
	    nodeHashList.add(fromColumnIndex + i, new HashMap<Integer, Node>);
	    Node cur = new Node('.', fromColumnIndex + i + 1);
	    this.addVertex(cur);//add vertex and add to nodeHashList
	    if(pre != null)
		this.g.addEdge(pre, cur);
	    gapNodes[i] = cur;
	    pre = cur;
	}
	*/
	//NEED TO SHIFT all columns after insertion, so updating all columnIndex (originalIndex+insSize.
	for(int i=fromColumnIndex+insSize; i<this.nodeHashList.size(); i++)
	    this.adjustColumnIndex(this.nodeHashList.get(i), i);//this.adjustColumnIndex(i);
	
	//remove all edges between start nodes and end nodes and add new edges connecting through gap nodes.
	double weightSum = this.getWeightSumsBetween2Columns(startNodes, endNodes, gapNodes);
	
	
	if(insSize > 1){
	    for(int i=fromColumnIndex; i<fromColumnIndex+insSize-1; i++){
		gapNodes = Arrays.copyOfRange(gapNodes, 1, gapNodes.length);
		this.getWeightSumsBetween2Columns(this.nodeHashList.get(i), endNodes, gapNodes);
	    }
	}
	
    }


    
    //removes all edges betweend start nodes and end nodes
    //connect edges to newly added gap nodes with correct weights
    private double getWeightSumsBetween2Columns(HashMap<Integer, Node> start,  HashMap<Integer, Node> end, Node[] gapNodes){
	Node sGap = gapNodes[0];
	Node eGap = gapNodes[gapNodes.length-1];
	
	double[] outweight = new double[6]; /* for each nucleotide */
	//ArrayList<Byte>[] outFScore = new ArrayList<Byte>[5];
	//ArrayList<Byte>[] outRScore = new ArrayList<Byte>[5];
	
	ArrayList<ArrayList<Byte>> outFScore = new ArrayList<ArrayList<Byte>>();
	ArrayList<ArrayList<Byte>> outRScore = new ArrayList<ArrayList<Byte>>();
	
	
	double[] inweight = new double[6];
	//ArrayList<Byte>[] inFScore = new ArrayList<Byte>[5];
	//ArrayList<Byte>[] inRScore = new ArrayList<Byte>[5];

	ArrayList<ArrayList<Byte>> inFScore = new ArrayList<ArrayList<Byte>>();
	ArrayList<ArrayList<Byte>> inRScore = new ArrayList<ArrayList<Byte>>();
	
	//ArrayList<HashSet<Integer>> outRHash = new ArrayList<HashSet<Integer>>();
	//ArrayList<HashSet<Integer>> inRHash = new ArrayList<HashSet<Integer>>();
	ArrayList<CustomHashMap> outRHash = new ArrayList<CustomHashMap>();
	ArrayList<CustomHashMap> inRHash = new ArrayList<CustomHashMap>();
	
	
	//for each nucleotide
	for(int i=0; i<6; i++){
	    outFScore.add(new ArrayList<Byte>());
	    outRScore.add(new ArrayList<Byte>());
	    inFScore.add(new ArrayList<Byte>());
	    inRScore.add(new ArrayList<Byte>());
	    //outRHash.add(new HashSet<Integer>());
	    //inRHash.add(new HashSet<Integer>());
	    outRHash.add(new CustomHashMap());
	    inRHash.add(new CustomHashMap());
	}
	
	double sum = 0.0d;
	//HashSet<Integer> rHashForGapNodes = new HashSet<Integer>();
	CustomHashMap rHashForGapNodes = new CustomHashMap();//new HashSet<Integer>();
	
	Integer[] sKeys = new Integer[0];
	Integer[] eKeys = new Integer[0];
	sKeys = start.keySet().toArray(sKeys);
	eKeys = end.keySet().toArray(eKeys);
  	/*
	for(int i=0;i<eKeys.length; i++){
	    rHashForGapNodes.addAll(end.get(eKeys[i].intValue()).getReadHashSet());
	    }*/
	
	boolean[] sEdgePresent = new boolean[6];
	boolean[] eEdgePresent = new boolean[6];
	boolean isThereConnection = false;
	
	//check all edges between starNodes and endNodes and sum up baseWise.
	for(int i=0; i < sKeys.length; i++){
	    int sVal = sKeys[i].intValue();
	    //if(sVal != 4){//edges between gap nodes are skipped, taken care of separately
		Node stNode = start.get(sKeys[i]);
		for(int j=0; j < eKeys.length; j++){
		    int eVal = eKeys[j].intValue();
		    //if(eVal != 4){//edges between gap nodes are skipped, taken care of separately
			Node eNode = end.get(eKeys[j]);
			CustomWeightedEdge e = this.g.getEdge(stNode, eNode);
			if(e != null){
			    sEdgePresent[sVal] = true;
			    eEdgePresent[eVal] = true;
			    isThereConnection = true;
			    double w = this.g.getEdgeWeight(e);
			    outweight[sVal] += w;
			    outFScore.get(sVal).addAll(e.getFScores());
			    outRScore.get(sVal).addAll(e.getRScores());
			    inweight[eVal] += w;
			    inFScore.get(eVal).addAll(e.getFScores());
			    inRScore.get(eVal).addAll(e.getRScores());
			    outRHash.get(sVal).addAll(e.getReadHashSet());
			    inRHash.get(eVal).addAll(e.getReadHashSet());
			    rHashForGapNodes.addAll(e.getReadHashSet());
			    sum += w;
			    this.g.removeEdge(e);
			}
			//		    }
		}
		//	    }
	}

	//we only need to add edges if there were edges between start and end
	if(isThereConnection){
	
	    //setting outgoing edges from start nodes to newly added gapNode( sGap ).
	    for(int i=0; i<sKeys.length; i++){
		if(sEdgePresent[sKeys[i].intValue()]){
		    Node stNode = start.get(sKeys[i]);
		    CustomWeightedEdge e = this.g.getEdge(stNode, sGap);
		    if(e == null){
			e = this.g.addEdge(stNode, sGap);
			this.g.setEdgeWeight(e, 0.0d);
		    }
		    this.g.setEdgeWeight(e, this.g.getEdgeWeight(e) + outweight[sKeys[i].intValue()]);//this.setEdgeWeight(e, outweight[sKeys[i].intValue()]);
		    e.addAllReadsFrom(outRHash.get(sKeys[i].intValue()));
		    e.addAllFScores(outFScore.get(sKeys[i].intValue()));
		    e.addAllRScores(outRScore.get(sKeys[i].intValue()));
		}
	    }
	    
	    //setting incoming edges from newly added gapNode( eGap ) to end nodes.
	    for(int i=0; i<eKeys.length; i++){
		if(eEdgePresent[eKeys[i].intValue()]){
		    Node eNode = end.get(eKeys[i]);
		    CustomWeightedEdge e = this.g.getEdge(eGap, eNode);//this.g.addEdge(eGap, eNode);
		    if(e == null){
			e = this.g.addEdge(eGap, eNode);
			this.g.setEdgeWeight(e, 0.0d);
		    }
		    this.g.setEdgeWeight(e, this.g.getEdgeWeight(e) + inweight[eKeys[i].intValue()]);
		    e.addAllReadsFrom(inRHash.get(eKeys[i].intValue()));
		    e.addAllFScores(inFScore.get(eKeys[i].intValue()));
		    e.addAllRScores(inRScore.get(eKeys[i].intValue()));
		}
	    }
	    
	    //set edgeWeight between newly inserted gap nodes.
	    //and add read identifiers to gapNodes
	    for(int i=0; i<gapNodes.length; i++){
		if(i>0){
		    CustomWeightedEdge e = this.g.getEdge(gapNodes[i-1], gapNodes[i]);
		    this.g.setEdgeWeight(e, this.g.getEdgeWeight(e) + sum);//this.g.getEdge(gapNodes[i-1], gapNodes[i]), sum);
		    e.addAllReadsFrom(rHashForGapNodes);
		}
		//gapNodes[i].addAllReadsFrom(rHashForGapNodes);
	    }
	}

	return sum;
    }

    //set columnIndex to newIndex.
    /*
    private void adjustColumnIndex(int newIndex){
	HashMap<Integer, Node> curHash = this.nodeHashList.get(newIndex);
	Iterator<Integer> keys = curHash.keySet().iterator();
	while(keys.hasNext())
	    curHash.get(keys.next()).setColIndex(newIndex);
    }
    */
    
    private void adjustColumnIndex(HashMap<Integer, Node> hash, int newIndex){
	Iterator<Integer> keys = hash.keySet().iterator();
	while(keys.hasNext())
	    hash.get(keys.next()).setColIndex(newIndex);
    }
    
    public void removeUnused(){
	this.removeUnusedEdges();
	this.removeUnusedVertices();
    }

    /* remove low frequency edges */
    private void removeUnusedEdges(){
	Iterator<CustomWeightedEdge> itr = this.g.edgeSet().iterator();
	CustomWeightedEdge e = null;
	ArrayList<CustomWeightedEdge> removalList = new ArrayList<CustomWeightedEdge>();
	
	while(itr.hasNext()){
	    e = itr.next();
	    if(this.g.getEdgeWeight(e) < 1.0d){
		removalList.add(e);//this.g.removeEdge(e);
	    }
	}
	System.err.println(this.HLAGeneName +"\t:removed\t" + removalList.size() + "\tEdges." );
	for(int i=0; i<removalList.size(); i++){
	    this.g.removeEdge(removalList.get(i));
	}
    }
    
    /* remove island vertices */
    private void removeUnusedVertices(){
	Iterator<Node> itr = this.g.vertexSet().iterator();
	Node n = null;
	ArrayList<Node> removalList = new ArrayList<Node>();
	while(itr.hasNext()){
	    n = itr.next();
	    //we dont remove sNode and tNode
	    if(!n.equals(this.sNode) && !n.equals(this.tNode)){
		if(this.g.inDegreeOf(n) ==0  && this.g.outDegreeOf(n) == 0){//this.g.degreeOf(n) < 1){
		    removalList.add(n);
		    //this.removeVertexFromNodeHashList(n);
		    //this.g.removeVertex(n);
		}
	    }
	}
	System.err.println(this.HLAGeneName +"\t:removed\t" + removalList.size() + "\tVertices." );
	for(int i=0; i<removalList.size(); i++){
	    //System.err.println("\t" + removalList.get(i).toString());
	    this.removeVertex(removalList.get(i));
	    //this.removeVertexFromNodeHashList(removalList.get(i));
	    //this.g.removeVertex(removalList.get(i));
	}
    }

    //removing stems. (unreachable stems and dead-end stems)
    /*
    public void removeStems(){
	Iterator<Node> itr= this.g.vertexSet().iterator();
	//Set<Node> vSet = this.g.vertexSet();
	//Node[] nodes = new Node[vSet.size()];
	Node n = null;
	ArrayList<Node> removalNodes = new ArrayList<Node>(); 
	
	int terminalStem = 0;
	int unreachableStem = 0;

	
	while(itr.hasNext()){
	    n = itr.next();
	    if(!n.equals(this.sNode) && !n.equals(this.tNode)){
		
		//dead-end stem    ---->x--->x
		if(this.g.outDegreeOf(n) == 0 && this.g.inDegreeOf(n) == 1){
		    int stemSize = 0;
		    terminalStem++;
		    Node curNode = n;
		    while(true){
			removalNodes.add(curNode);
			stemSize++;
			CustomWeightedEdge e = this.g.incomingEdgesOf(curNode).toArray(new CustomWeightedEdge[1])[0];
			System.err.print(this.g.getEdgeWeight(e) + "\t");
			Node nextNode = this.g.getEdgeSource(e);
			if(this.g.outDegreeOf(nextNode) == 1 && this.g.inDegreeOf(nextNode) == 1)
			    curNode = nextNode;
			else
			    break;
			
		    }
		    System.err.println();
		    System.err.println("[DE]stemSize:\t" + stemSize);
		}
		//unreachable stem   x--->x--->
		else if(this.g.outDegreeOf(n) == 1 && this.g.inDegreeOf(n) == 0){
		    int stemSize = 0;
		    unreachableStem++;
		    Node curNode = n;
		    while(true){
			removalNodes.add(curNode);
			stemSize++;
			CustomWeightedEdge e = this.g.outgoingEdgesOf(curNode).toArray(new CustomWeightedEdge[1])[0];
			System.err.print(this.g.getEdgeWeight(e) + "\t");
			Node nextNode = this.g.getEdgeSource(e);
			if(this.g.outDegreeOf(nextNode) == 1 && this.g.inDegreeOf(nextNode) == 1)
			    curNode = nextNode;
			else
			    break;
			
		    }
		    System.err.println("[UN]stemSize:\t" + stemSize);
		}
	    }
	}
	System.err.println(this.HLAGeneName + "\t:removed\t[DE]:" + terminalStem + "\t[UN]:" + unreachableStem + "\t[NumVertices]:" + removalNodes.size());
	for(int i=0; i<removalNodes.size(); i++){
	    this.removeVertex(removalNodes.get(i));
	    //this.removeVertexFromNodeHashList(removalNodes.get(i));
	    //this.g.removeVertex(removalNodes.get(i));
	}
    }
    */
    
    /* remove any stems */
    public void removeStems(){
	ArrayList<int[]> typingIntervals = this.obtainTypingIntervals();
	
	//Set<Node> vSet = this.g.vertexSet();
	Node[] nodes = this.g.vertexSet().toArray(new Node[0]);//new Node[vSet.size()];
	HashSet<Node> dNodes = new HashSet<Node>();
	Node n = null;
	int terminalStem = 0;
	int unreachableStem = 0;
	for(int i=0; i<nodes.length; i++){
	    n = nodes[i];
	    if(!n.equals(this.sNode) && !n.equals(this.tNode) && !dNodes.contains(n)){
		
		//dead-end stem    ---->x--->x
		if(this.g.outDegreeOf(n) == 0 && this.g.inDegreeOf(n) == 1){
		    int stemSize = 0;
		    terminalStem++;
		    Node curNode = n;
		    
		    while(true){
			if(!this.alleles.get(0).withinTypingRegion(curNode, typingIntervals))
			    ;//System.err.println("NOT IN TYPING INTERVAL!!");
			else
			    System.err.print("YES! IN TYPING INTERVAL!!");
			stemSize++;
			CustomWeightedEdge e = this.g.incomingEdgesOf(curNode).toArray(new CustomWeightedEdge[1])[0];
			System.err.print("\t" + this.g.getEdgeWeight(e));
			Node nextNode = this.g.getEdgeSource(e);
			dNodes.add(curNode);
			this.removeVertex(curNode);
			if(this.g.outDegreeOf(nextNode) == 0 && this.g.inDegreeOf(nextNode) == 1)
			    curNode = nextNode;
			else
			    break;
		    }
		    System.err.println("[DE]stemSize:\t" + stemSize);
		}
		//unreachable stem   x--->x--->
		else if(this.g.outDegreeOf(n) == 1 && this.g.inDegreeOf(n) == 0){
		    int stemSize = 0;
		    unreachableStem++;
		    Node curNode = n;
		    while(true){
			if(!this.alleles.get(0).withinTypingRegion(curNode, typingIntervals))
			    ;//System.err.println("NOT IN TYPING INTERVAL!!");
			else
			    System.err.println("YES! IN TYPING INTERVAL!!");
			stemSize++;
			CustomWeightedEdge e = this.g.outgoingEdgesOf(curNode).toArray(new CustomWeightedEdge[1])[0];
			System.err.print("\t" + this.g.getEdgeWeight(e));
			Node nextNode = this.g.getEdgeTarget(e);
			dNodes.add(curNode);
			this.removeVertex(curNode);
			if(this.g.outDegreeOf(nextNode) == 1 && this.g.inDegreeOf(nextNode) == 0)
			    curNode = nextNode;
			else
			    break;
		    }
		    System.err.println("[UN]stemSize:\t" + stemSize);
		}
	    }
	}
	System.err.println(this.HLAGeneName + "\t:removed\t[DE]:" + terminalStem + "\t[UN]:" + unreachableStem + "\t[NumVertices]:" + dNodes.size());
    }

    

    public void countStems(){
    
	Iterator<Node> itr = this.g.vertexSet().iterator();
	Node n = null;
	int terminalType = 0;
	int startType = 0;
	while(itr.hasNext()){
	    n = itr.next();
	    if(!n.equals(this.sNode) && !n.equals(this.tNode)){
		if(this.g.inDegreeOf(n) == 1 && this.g.outDegreeOf(n) == 0){
		    terminalType++;
		}else if(this.g.inDegreeOf(n) == 0 && this.g.outDegreeOf(n) == 1){
		    startType++;
		    System.err.println("startType:\t" + n.toString());
		}
	    }
	}
	System.err.println("Stems\t" + terminalType + "\t" + startType);
    }
    
    
    /*
    private void initNumPathForColumn(HashMap){
    
    }*/
    /*
    private ArrayList<Integer> getNextEdgeEncodingNumbersAtColumnC(int c){
	HashMap<Integer, Node> h = this.nodeHashList.get(c);
	Set<Integer> s = h.keySet();
	
	}*/


    
}


