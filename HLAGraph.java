import htsjdk.samtools.SAMRecord;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Set;
import java.util.Arrays;
import java.util.HashSet;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;

import org.jgrapht.*;
import org.jgrapht.graph.*;

public class HLAGraph{

    private String HLAGeneName; //A B C ... //HLA-DPA1, HLA-DPB1, HLA-DQA1, HLA-DQB1, HLA-DRA, and HLA-DRB1
    private ArrayList<Sequence> alleles; //
    private HashMap<String, Sequence> alleleHash;
    
    //private SimpleDirectedWeightedGraph<Node, DefaultWeightedEdge> g;
    private SimpleDirectedWeightedGraph<Node, CustomWeightedEdge> g;

    //private ArrayList<HashMap<Character, Node>> nodeHashList;//list index = columnIndex-1.
    private ArrayList<HashMap<Integer, Node>> nodeHashList;// list index = columnIndex - 1;

    private Node sNode;
    private Node tNode;
    
    private int columnLen;
    
    /* Outer list index = columnIndex -1 --> insertion point */
    /* Inner list index insertion length */ 
    //private ArrayList<ArrayList<HashMap<Character, Node>>> insertionNodeHashList;
    private ArrayList<ArrayList<HashMap<Integer, Node>>> insertionNodeHashList;


    /* DO NOT use this to add sNode and tNode */
    public void addVertex(Node n){
	int ci = n.getColIndex();
	this.g.addVertex(n);
	this.nodeHashList.get(n.getColIndex()-1).put(new Integer(n.getIBase()), n);
    }
    
    /* DO NOT use this to add sNode and tNode */
    public void removeVertex(Node n){
	this.nodeHashList.get(n.getColIndex()-1).remove(new Integer(n.getIBase()));
	this.g.removeVertex(n);
    }


    public void setHLAGeneName(String gn){
	this.HLAGeneName = gn;
    }
    
    
    public Sequence getRefAllele(){
	return this.alleles.get(0);
    }

    public HLAGraph(ArrayList<Sequence> seqs){
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

    
    //modified so that if pre node is null, create curnode but dont' attempt to connect w/ an edge
    private Node addMissingNode(char b, int colPos, Node cur, Node pre, boolean isRefStrand, byte qual, int readNum){
	cur = new Node(b, colPos);
	//cur.addRead(readNum);
	this.g.addVertex(cur);
	//this.nodeHashList.get(colPos - 1).put(new Character(b), cur);
	this.nodeHashList.get(colPos - 1).put(new Integer(Base.char2ibase(b)), cur);
	if(pre != null){
	    //DefaultWeightedEdge e = this.g.addEdge(pre, cur);
	    this.addAndIncrement(pre, cur, isRefStrand, qual, readNum);
	}else
	    cur.addRead(readNum);
	return cur;
    }

    private void addAndIncrement(Node source, Node target, boolean isRefStrand, byte qual, int readNum){
	target.addRead(readNum);
	CustomWeightedEdge e = this.g.addEdge(source,target);
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
	    target.addRead(readNum);
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
    

    public boolean traverseAndWeights(){
	System.err.println("=========================");
	System.err.println("=  " + this.HLAGeneName);
	System.err.println("=========================");


	Sequence ref = this.alleles.get(0);
	int fCount = 0;
	/* temporary typing interval information */
	ArrayList<int[]> typingIntervals = new ArrayList<int[]>();
	if(this.isClassI()){
	    int[] tmp = new int[2];
	    tmp[0] = ref.getBoundaries()[3];
	    tmp[1] = ref.getBoundaries()[4];
	    
	    typingIntervals.add(tmp);
	    
	    int[] tmp2 = new int[2];
	    tmp2[0] = ref.getBoundaries()[5];
	    tmp2[1] = ref.getBoundaries()[6];
	    
	    typingIntervals.add(tmp2);
	}else if(this.isClassII()){
	    int[] tmp2 = new int[2];
	    tmp2[0] = ref.getBoundaries()[3];
	    tmp2[1] = ref.getBoundaries()[4];
	    
	    typingIntervals.add(tmp2);
	}
	
	
	
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
		System.err.println("\nNEXT_EXON\n");
		for(int k=start-2; k<end; k++){
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
	    System.err.println(e.toString());
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
	
	Sequence ref = this.alleles.get(0);
	int fCount = 0;
	/* temporary typing interval information */
	ArrayList<int[]> typingIntervals = new ArrayList<int[]>();
	if(this.isClassI()){
	    int[] tmp = new int[2];
	    tmp[0] = ref.getBoundaries()[5];
	    tmp[1] = ref.getBoundaries()[6];
	    
	    typingIntervals.add(tmp);
	    
	    int[] tmp2 = new int[2];
	    tmp2[0] = ref.getBoundaries()[3];
	    tmp2[1] = ref.getBoundaries()[4];
	    
	    typingIntervals.add(tmp2);
	}else if(this.isClassII()){
	    int[] tmp2 = new int[2];
	    tmp2[0] = ref.getBoundaries()[3];
	    tmp2[1] = ref.getBoundaries()[4];
	    
	    typingIntervals.add(tmp2);
	}
	
	
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
	Integer[] sKeys = new Integer[5];
	sKeys = s.keySet().toArray(sKeys);
	Integer[] eKeys = new Integer[5];
	eKeys = s.keySet().toArray(eKeys);
	
	for(int i=0;i<sKeys.length; i++){
	    if(sKeys[i].intValue() != 4){
		for(int j=0; j<eKeys.length; j++){
		    if(eKeys[j].intValue() != 4){
			CustomWeightedEdge e = this.g.getEdge(s.get(sKeys[i]), t.get(eKeys[i]));
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
	
	for(int i=0; i<5; i++){
	    outFScore.add(new ArrayList<Byte>());
	    outRScore.add(new ArrayList<Byte>());
	    inFScore.add(new ArrayList<Byte>());
	    inRScore.add(new ArrayList<Byte>());
	}
	
	double sum = 0.0d;
	HashSet<Integer> rHashForGapNodes = new HashSet<Integer>();
	
	Integer[] sKeys = new Integer[5];
	Integer[] eKeys = new Integer[5];
	sKeys = start.keySet().toArray(sKeys);
	eKeys = end.keySet().toArray(eKeys);
  	
	for(int i=0;i<eKeys.length; i++){
	    rHashForGapNodes.addAll(end.get(eKeys[i].intValue()).getReadHashSet());
	}
	
	boolean[] sEdgePresent = new boolean[5];
	boolean[] eEdgePresent = new boolean[5];
	boolean isThereConnection = false;
	
	//check all edges between starNodes and endNodes and sum up baseWise.
	for(int i=0; i < sKeys.length; i++){
	    int sVal = sKeys[i].intValue();
	    if(sVal != 4){//edges between gap nodes are skipped, taken care of separately
		Node sNode = start.get(sKeys[i]);
		for(int j=0; j < eKeys.length; j++){
		    int eVal = eKeys[j].intValue();
		    if(eVal != 4){//edges between gap nodes are skipped, taken care of separately
			Node eNode = end.get(eKeys[j]);
			CustomWeightedEdge e = this.g.getEdge(sNode, eNode);
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
			    sum += w;
			    this.g.removeEdge(e);
			}
		    }
		}
	    }
	}

	//we only need to add edges if there were no edges between start and end
	if(isThereConnection){
	
	    //setting outgoing edges from start nodes to newly added gapNode( sGap ).
	    for(int i=0; i<sKeys.length; i++){
		if(sEdgePresent[sKeys[i].intValue()]){
		    Node sNode = start.get(sKeys[i]);
		    CustomWeightedEdge e = this.g.getEdge(sNode, sGap);
		    if(e == null){
			e = this.g.addEdge(sNode, sGap);
			this.g.setEdgeWeight(e, 0.0d);
		    }
		    this.g.setEdgeWeight(e, this.g.getEdgeWeight(e) + outweight[sKeys[i].intValue()]);//this.setEdgeWeight(e, outweight[sKeys[i].intValue()]);
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
		}
		gapNodes[i].addAllReadsFrom(rHashForGapNodes);
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

    private void removeUnusedEdges(){
	Iterator<CustomWeightedEdge> itr = this.g.edgeSet().iterator();
	CustomWeightedEdge e = null;
	ArrayList<CustomWeightedEdge> removalList = new ArrayList<CustomWeightedEdge>();
	while(itr.hasNext()){
	    e = itr.next();
	    if(this.g.getEdgeWeight(e) < 2.0d){
		removalList.add(e);//this.g.removeEdge(e);
	    }
	}
	System.err.println(this.HLAGeneName +"\t:removed\t" + removalList.size() + "\tEdges." );
	for(int i=0; i<removalList.size(); i++){
	    this.g.removeEdge(removalList.get(i));
	}
    }
    
    private void removeUnusedVertices(){
	Iterator<Node> itr = this.g.vertexSet().iterator();
	Node n = null;
	ArrayList<Node> removalList = new ArrayList<Node>();
	while(itr.hasNext()){
	    n = itr.next();
	    //we dont remove sNode and tNode
	    if(!n.equals(this.sNode) && !n.equals(this.tNode)){
		if(this.g.inDegreeOf(n) < 1 || this.g.outDegreeOf(n) < 1){//this.g.degreeOf(n) < 1){
		    removalList.add(n);
		    //this.removeVertexFromNodeHashList(n);
		    //this.g.removeVertex(n);
		}
	    }
	}
	System.err.println(this.HLAGeneName +"\t:removed\t" + removalList.size() + "\tVertices." );
	for(int i=0; i<removalList.size(); i++){
	    System.err.println("\t" + removalList.get(i).toString());
	    this.removeVertexFromNodeHashList(removalList.get(i));
	    this.g.removeVertex(removalList.get(i));
	}
    }
    
    //removes node from nodeHashList. We dont touch insertionNodeHashList because any node added on insertionNodeHashList must have weights.
    private void removeVertexFromNodeHashList(Node n){
	this.nodeHashList.get(n.getColIndex()-1).remove(new Integer(n.getIBase()));
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


