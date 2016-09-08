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
    private ArrayList<Path> interBubblePaths;

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
    /*
    //first exon of typing region
    private int headerExcessLengthBeyondTypingBoundary1; 
    private int tailExcessLengthBeyondTypingBoundary1;
    
    //secon exon of typing region
    private int headerExcessLengthBeyondTypingBoundary2;
    private int tailExcessLengthBeyondTypingBoundary2;
    */  
    /* Outer list index = columnIndex -1 --> insertion point */
    /* Inner list index insertion length */ 
    //private ArrayList<ArrayList<HashMap<Character, Node>>> insertionNodeHashList;
    private ArrayList<ArrayList<HashMap<Integer, Node>>> insertionNodeHashList;

    public void setTypingSequences(ArrayList<HLASequence> seqs){
	this.typingSequences = seqs;
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
	//this.headerExcessLengthBeyondTypingBoundary2 = 0;
	//this.tailExcessLengthBeyondTypingBoundary2 = 0;
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
	e.addRead(readNum);
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
	    e.addRead(readNum);
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
			//colPos = curAllele.getColPosFromBasePos(refBasePos);
			for(int i=0; i<cigarLen; i++){
			    numOp++;
			    int tmpColPos = curAllele.getNextColPosForBase(colPos - 1) + 1;
			    if(tmpColPos > colPos){
				for(int j=colPos;j<tmpColPos;j++){
				    HLA.HOPPING++;
				    curnode = this.nodeHashList.get(j-1).get(new Integer(Base.char2ibase('.')));
				    //this.incrementWeight(prevnode,curnode);
				    this.incrementWeight(prevnode,curnode,isRefStrand, quals[baseIndex-1], readNum);
				    prevnode=curnode;
				}
				colPos = tmpColPos;
			    }

			    //colPos = curAllele.getNextcolPosForBase(colpos - 1);
			    //System.err.println("Processing position : " + i + "(colPos:" + colPos + "[" + (char)bases[baseIndex] + "])");//+ Base.ibase2char(bases[baseIndex]) + "])");
			    //int colPos = curAllele.getColPosFromBasePos(refBasePos);
			    //System.out.println((curnode==null ? "curnode<null>" : "curnode<NOTnull>") + "\t" + (prevnode==null ? "prevnode<null>" : "prevnode<NOTnull>"));
			    //curnode = this.nodeHashList.get(colPos -1).get(new Character((char)bases[baseIndex]));
			    curnode = this.nodeHashList.get(colPos -1).get(new Integer(Base.char2ibase((char)bases[baseIndex])));
			    
			    //System.err.println((curnode==null ? "curnode<null>" : "curnode<NOTnull>") + "\t" + (prevnode==null ? "prevnode<null>" : "prevnode<NOTnull>"));
			    //if no such node found, we add new node and add edge from prevnode.
			    if(curnode == null){
				

				// THIS CASE IS HANDLED BY modifying addMissingNode.
				//this means it starts with mismatch 
				//if(prevnode == null)
				//  System.err.println("Handle This?????");
				
				/*
				Set<Integer> keys = this.nodeHashList.get(colPos - 1).keySet();
				Iterator<Integer> itr = keys.iterator();
				System.err.print("colPos[" + colPos + "]:\t");
				while(itr.hasNext()){
				    Integer tmpI = itr.next();
				    System.err.print(Base.ibase2char(tmpI.intValue()) + "\t");
				}
				System.err.println();
				*/
				HLA.NEW_NODE_ADDED++;
				//curnode = this.addMissingNode((char)bases[baseIndex], colPos, curnode, prevnode);
				curnode = this.addMissingNode((char)bases[baseIndex], colPos, curnode, prevnode, isRefStrand, quals[baseIndex], readNum);
				if(curnode == null)
				    System.err.println("IMPOSSIBLE: curnode NULL again after adding missing node!");
			    }
			    else if(prevnode != null){
				//this.incrementWeight(prevnode, curnode);//source, target);
				this.incrementWeight(prevnode, curnode, isRefStrand, quals[baseIndex], readNum);
			    }
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
			    int tmpColPos = curAllele.getNextColPosForBase(colPos - 1) + 1;
			    if(tmpColPos > colPos){
				for(int j=colPos;j<tmpColPos;j++){
				    HLA.HOPPING++;
				    curnode = this.nodeHashList.get(j-1).get(new Integer(Base.char2ibase('.')));
				    //this.incrementWeight(prevnode,curnode);
				    this.incrementWeight(prevnode, curnode, isRefStrand, quals[baseIndex-1], readNum);
				    prevnode=curnode;
				}
				colPos = tmpColPos;
			    }
			    
			    //curnode = this.nodeHashList.get(colPos - 1).get(new Character('.'));
			    curnode = this.nodeHashList.get(colPos - 1).get(new Integer(Base.char2ibase('.')));
			    //System.err.println((curnode==null ? "curnode<null>" : "curnode<NOTnull>") + "\t" + (prevnode==null ? "prevnode<null>" : "prevnode<NOTnull>"));
			    if(curnode == null){
				HLA.NEW_NODE_ADDED++;
				//System.err.println("HERE (D)");
				//curnode = this.addMissingNode('.', colPos, curnode, prevnode);
				curnode = this.addMissingNode('.', colPos, curnode, prevnode, isRefStrand, quals[baseIndex-1], readNum);
			    }else{
				//this.incrementWeight(prevnode, curnode);
				this.incrementWeight(prevnode, curnode, isRefStrand, quals[baseIndex-1], readNum);
			    }
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
				    System.err.println("SHOULD NOT HAPPEND (2)[addWeight]");
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
	for(int i=0;i<5;i++){
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
	for(int i=0; i<bubbles.size(); i++){
	    bubbles.get(i).initBubbleSequences();
	}
	
	/* superBubble is a merged bubbles. Ideally, you want to have just one bubble. */
	ArrayList<Bubble> superBubbles = new ArrayList<Bubble>();
	
	
	Bubble curSuperBubble = bubbles.get(0);
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
	    MergeStatus ms = curSuperBubble.mergeBubble(bubbles.get(i), lastSegregationColumnIndex);
	    
	    //if we are cutting here
	    if(ms.isSplit()){
		System.out.println("CANT PHASE --> setting OB as curSuperBubble.");
		superBubbles.add(curSuperBubble);
		curSuperBubble = bubbles.get(i);
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
	
	this.printBubbleResults(superBubbles);
	this.getFracturedPaths(superBubbles, this.headerExcessLengthBeyondTypingBoundary, this.tailExcessLengthBeyondTypingBoundary);
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
    public void printBubbleResults(ArrayList<Bubble> superBubbles){
	//StringBuffer output = new StringBuffer();
	int startIndex = 0;
	
	System.out.println("Printing\t" + superBubbles.size() + "\tfractured super bubbles.");
	//output.append(superBubbles.size() + "\tfractured SuperBubbles\n");
	int count = 0;


	ArrayList<ArrayList<DNAString>> fracturedSequences = new ArrayList<ArrayList<DNAString>>();

	for(Bubble sb : superBubbles){
	    ArrayList<DNAString> sequences = new ArrayList<DNAString>();
	    fracturedSequences.add(sequences);
	    System.out.println("\tSuperBubble\t" + count);
	    startIndex = sb.printResults(this.interBubbleSequences, startIndex, sequences, this.HLAGeneName , count);
	    count++;
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

	this.generateCandidates(fracturedSequences);
    }


    public void getFracturedPaths(ArrayList<Bubble> superBubbles, int[] headerExcessArr, int[] tailExcessArr){
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
		/*if(subj.getGroup().getGroupString().equals("A*01:01:01G")){
		    System.err.println(candidate);
		    System.err.println(subject);
		    System.err.println("A*01:01:01G\t" + curR.toString());
		    }*/
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
	ArrayList<Integer> bubbleLengths = new ArrayList<Integer>();
	ArrayList<Integer> coordinates = new ArrayList<Integer>(); //keeps track of start coordinates of bubbles
	/* counters */
	
	Node curSNode = null;

	this.interBubbleSequences = new ArrayList<StringBuffer>();
	this.interBubblePaths = new ArrayList<Path>();
	
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

	    /* FOR EACH POSITION in a TYPING INTERVAL*/
	    for(int k=start-1;k<end-1;k++){
		HashMap<Integer, Node> columnHash = this.nodeHashList.get(k);
		Integer[] keys = columnHash.keySet().toArray(new Integer[0]);
		
		/*it's a collapsing node if curBubbleLength > 2
		  else it's a possible start of bubble.*/
		if(keys.length == 1){
		    //headerBubble = false;
		    /* then it must be a collapsing node; */
		    if(curBubbleLength > 1){
			this.interBubbleSequences.add(curbf);
			this.interBubblePaths.add(tp.toPath(this.g));
			//this.interBubblePaths.add(curP);
			curBubbleLength++;
			numBubbles++;
			//numPaths.add(new Integer(this.analyzeBubble(lastStartOfBubble, k)));
			bubbleLengths.add(new Integer(curBubbleLength-2));
			coordinates.add(new Integer(lastStartOfBubble));
			if(firstBubble){
			    bubbles.add(new Bubble(this, curSNode, columnHash.get(keys[0]), firstBubble));
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
		}else if(keys.length > 1){

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
				curbf.append(curSNode.getBase());
				tp.appendNode(curSNode);
				lastStartOfBubble = l;
				curBubbleLength = tmpBubbleLength;
				this.headerExcessLengthBeyondTypingBoundary[i] = curBubbleLength - 1;
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
		    }else{
			curBubbleLength++;
			//preNode = null;
		    }
		}else{//disconnected graph.
		    System.err.println("This should NOT HAPPEN");
		}
	    }
	    this.interBubbleSequences.add(curbf);
	    this.interBubblePaths.add(tp.toPath(this.g));
	    curbf = new StringBuffer("");
	    tp = new TmpPath();
	    if(curBubbleLength > 1){
		System.err.println(">>>>>>>Bubble at the end:\t[curBubbleLength]:"+ curBubbleLength);
	    }
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


    public void flattenInsertionNodes(){
	
	ArrayList<int[]> typingIntervals = this.obtainTypingIntervals();

	int fCount = 0;	
	
	for(int i=typingIntervals.size()-1; i>-1; i--){
	    int start = typingIntervals.get(i)[0];
	    int end   = typingIntervals.get(i)[1];
	    
	    for(int j=end-1; j >= start; j--){
		int insSize = this.insertionNodeHashList.get(j).size();
		//there is insertion, we need to flatten.
		if(insSize > 0 && this.isThereConnectionToInsertionNodes(insSize, j+1)){
		    fCount++;
		    this.shiftColumnsByInsertionSize(insSize, j+1);
		}
	    }
	}
	
	System.err.println(this.HLAGeneName + "\t>>>>> FLATTENED InsertionBubble:\t" + fCount );
    }

    
    private boolean isThereConnectionToInsertionNodes(int insSize, int fromColumnIndex){
	HashMap<Integer, Node> startNodes = nodeHashList.get(fromColumnIndex-1);
	boolean sConnection = false;
	boolean eConnection = false;
	HashMap<Integer, Node> sInsHash = this.insertionNodeHashList.get(fromColumnIndex - 1).get(0);
	HashMap<Integer, Node> eInsHash = this.insertionNodeHashList.get(fromColumnIndex - 1).get(insSize - 1);
	HashMap<Integer, Node> endNodes = nodeHashList.get(fromColumnIndex);
	
	sConnection = this.isThereConnection(startNodes, sInsHash);
	eConnection = this.isThereConnection(eInsHash, endNodes);
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
		    System.err.print("eKyes[j] intval\t");
		    System.err.println(eKeys[j].intValue());
		    if(eKeys[j].intValue() != 4){
			CustomWeightedEdge e = this.g.getEdge(s.get(sKeys[i]), t.get(eKeys[j]));
			if(e != null)
			    return true;
		    }
		}
	    }		
	}
	return false;
    }

    /* fromColumnIndex is 0-based index */
    /* 0based(List index): 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 */
    /* 1based(CI in Node): 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 */
    /* from ColumnIndex at 5, insSize of 2*/
    
    private void shiftColumnsByInsertionSize(int insSize, int fromColumnIndex){
	
	HashMap<Integer, Node> startNodes = nodeHashList.get(fromColumnIndex-1);
	HashMap<Integer, Node> endNodes = nodeHashList.get(fromColumnIndex);

	//we need to insert <insSize>-many columns first
	Node pre = null;
	Node[] gapNodes = new Node[insSize];
	ArrayList<Base> insBases = new ArrayList<Base>();//new Base[insSize];


	//insert insSize-many columns with gapNodes and transfer insertionNodes to nodeHashList.
	for(int i=0; i<insSize; i++){
	    //add a space first then add the vertex --> gets the space(HashMap) from insertionNodeHashList
	    HashMap<Integer, Node> insHash_i = this.insertionNodeHashList.get(fromColumnIndex-1).get(i);
	    this.adjustColumnIndex(insHash_i, fromColumnIndex + i + 1);
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
	    this.alleles.get(i).insertBlanks(fromColumnIndex, insBases);
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
	
	double[] outweight = new double[5];
	//ArrayList<Byte>[] outFScore = new ArrayList<Byte>[5];
	//ArrayList<Byte>[] outRScore = new ArrayList<Byte>[5];
	
	ArrayList<ArrayList<Byte>> outFScore = new ArrayList<ArrayList<Byte>>();
	ArrayList<ArrayList<Byte>> outRScore = new ArrayList<ArrayList<Byte>>();
	
	
	double[] inweight = new double[5];
	//ArrayList<Byte>[] inFScore = new ArrayList<Byte>[5];
	//ArrayList<Byte>[] inRScore = new ArrayList<Byte>[5];

	ArrayList<ArrayList<Byte>> inFScore = new ArrayList<ArrayList<Byte>>();
	ArrayList<ArrayList<Byte>> inRScore = new ArrayList<ArrayList<Byte>>();
	
	ArrayList<HashSet<Integer>> outRHash = new ArrayList<HashSet<Integer>>();
	ArrayList<HashSet<Integer>> inRHash = new ArrayList<HashSet<Integer>>();
	
	for(int i=0; i<5; i++){
	    outFScore.add(new ArrayList<Byte>());
	    outRScore.add(new ArrayList<Byte>());
	    inFScore.add(new ArrayList<Byte>());
	    inRScore.add(new ArrayList<Byte>());
	    outRHash.add(new HashSet<Integer>());
	    inRHash.add(new HashSet<Integer>());
	}
	
	double sum = 0.0d;
	HashSet<Integer> rHashForGapNodes = new HashSet<Integer>();
	
	Integer[] sKeys = new Integer[0];
	Integer[] eKeys = new Integer[0];
	sKeys = start.keySet().toArray(sKeys);
	eKeys = end.keySet().toArray(eKeys);
  	/*
	for(int i=0;i<eKeys.length; i++){
	    rHashForGapNodes.addAll(end.get(eKeys[i].intValue()).getReadHashSet());
	    }*/
	
	boolean[] sEdgePresent = new boolean[5];
	boolean[] eEdgePresent = new boolean[5];
	boolean isThereConnection = false;
	
	//check all edges between starNodes and endNodes and sum up baseWise.
	for(int i=0; i < sKeys.length; i++){
	    int sVal = sKeys[i].intValue();
	    if(sVal != 4){//edges between gap nodes are skipped, taken care of separately
		Node stNode = start.get(sKeys[i]);
		for(int j=0; j < eKeys.length; j++){
		    int eVal = eKeys[j].intValue();
		    if(eVal != 4){//edges between gap nodes are skipped, taken care of separately
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
		    }
		}
	    }
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


