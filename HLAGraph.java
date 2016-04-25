import htsjdk.samtools.SAMRecord;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Set;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;

import org.jgrapht.*;
import org.jgrapht.graph.*;

public class HLAGraph{

    private String HLAGeneName; //A B C ...
    private ArrayList<Sequence> alleles; //
    private HashMap<String, Sequence> alleleHash;
    
    private SimpleDirectedWeightedGraph<Node, DefaultWeightedEdge> g;

    //private ArrayList<HashMap<Character, Node>> nodeHashList;//list index = columnIndex-1.
    private ArrayList<HashMap<Integer, Node>> nodeHashList;// list index = columnIndex - 1;

    private Node sNode;
    private Node tNode;
    
    /* Outer list index = columnIndex -1 --> insertion point */
    /* Inner list index insertion length */ 
    //private ArrayList<ArrayList<HashMap<Character, Node>>> insertionNodeHashList;
    private ArrayList<ArrayList<HashMap<Integer, Node>>> insertionNodeHashList;
    
    
    public HLAGraph(ArrayList<Sequence> seqs){
	this.alleles = seqs; 
	this.alleleHash = new HashMap<String, Sequence>();
	for(int i=0;i<this.alleles.size();i++){
	    this.alleleHash.put(this.alleles.get(i).getAlleleName(), this.alleles.get(i));
	}
	this.g = new SimpleDirectedWeightedGraph<Node, DefaultWeightedEdge>(DefaultWeightedEdge.class);
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
    private Node addMissingNode(char b, int colPos, Node cur, Node pre){
	cur = new Node(b, colPos);
	this.g.addVertex(cur);
	//this.nodeHashList.get(colPos - 1).put(new Character(b), cur);
	this.nodeHashList.get(colPos - 1).put(new Integer(Base.char2ibase(b)), cur);
	if(pre != null){
	    DefaultWeightedEdge e = this.g.addEdge(pre, cur);
	    this.g.setEdgeWeight(e, 1.0d);
	}
	return cur;
    }


    private void incrementWeight(Node source, Node target){
	DefaultWeightedEdge e = g.getEdge(source, target);
	if(e == null){
	    e = g.addEdge(source, target);
	    g.setEdgeWeight(e, 1.0d);
	}else
	    g.setEdgeWeight(e, g.getEdgeWeight(e)+1);
    }
    
    
    public int addWeight(SAMRecord sr){
	int numOp = 0;
	Cigar cigar = sr.getCigar();
	byte[] bases = sr.getReadBases(); //ASCII bytes ACGTN=.
	int baseIndex = 0;
	int refBasePos = sr.getAlignmentStart();
	Node prevnode = null;
	Node curnode = null;
	Base curbase = null;
	Sequence curAllele = this.alleleHash.get(sr.getReferenceName());
	int colPos = curAllele.getColPosFromBasePos(refBasePos);

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
				    this.incrementWeight(prevnode,curnode);
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
				curnode = this.addMissingNode((char)bases[baseIndex], colPos, curnode, prevnode);
				if(curnode == null)
				    System.err.println("IMPOSSIBLE: curnode NULL again after adding missing node!");
			    }
			    else if(prevnode != null)
				this.incrementWeight(prevnode, curnode);//source, target);
			    
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
				    this.incrementWeight(prevnode,curnode);
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
				curnode = this.addMissingNode('.', colPos, curnode, prevnode);
			    }else
				this.incrementWeight(prevnode, curnode);
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
				    
				    DefaultWeightedEdge e = this.g.addEdge(prevnode, curnode);
				    this.g.setEdgeWeight(e, 1.0d);
				}else
				    this.incrementWeight(prevnode, curnode);
				prevnode = curnode;
				baseIndex++;
			    }else if(tmpColPos > colPos){//then we must insert here.
				curnode = this.nodeHashList.get(colPos - 1).get(new Integer(Base.char2ibase((char)bases[baseIndex])));
				if(curnode == null){
				    HLA.NEW_NODE_ADDED++;
				    curnode = this.addMissingNode((char)bases[baseIndex], colPos, curnode, prevnode);
				    if(curnode == null){
					System.err.println("IMPOSSIBLE: curnode NULL again after adding missing node! (1)[addWeight]");
					System.exit(9);
				    }
				}else if(prevnode !=null){
				    HLA.INSERTION_WITH_NO_NEW_NODE++;
				    this.incrementWeight(prevnode, curnode);
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
		DefaultWeightedEdge e;
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

    
    public void traverseAndWeights(){
	Node preNode;
	Node curNode;
	for(int i=0; i<this.alleles.size(); i++){
	    preNode = this.sNode;
	    Sequence curseq = this.alleles.get(i);
	    double sum = 0;
	    int numZero = 0;
	    System.err.println(curseq.getAlleleName());
	    for(int j=0; j<curseq.getColLength(); j++){
		char uchar = Character.toUpperCase(curseq.baseAt(j).getBase());
		HashMap<Integer, Node> curHash = this.nodeHashList.get(j);
		curNode = this.nodeHashList.get(j).get(new Integer(Base.char2ibase(uchar)));
		if(!preNode.equals(this.sNode)){
		    System.err.print(uchar + "[" + this.g.getEdgeWeight(this.g.getEdge(preNode, curNode)) + "]->");
		    double tmpw = this.g.getEdgeWeight(this.g.getEdge(preNode, curNode));
		    sum+=tmpw;
		    if(tmpw == 0.0d)
			numZero++;
		    
		}
		preNode = curNode;
	    }
	    System.err.println("\n" + curseq.getAlleleName() + "\tSUM:\t" + sum + "\t#ZERO:\t" + numZero);
	}
    }

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

    
}


