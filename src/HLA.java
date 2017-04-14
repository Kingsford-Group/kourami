import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

import java.io.File;
import java.io.IOException;
import java.io.BufferedWriter;
import java.io.FileWriter;
import java.util.HashMap;
import java.util.Iterator;
import java.util.ArrayList;

import it.unimi.dsi.fastutil.objects.Object2IntOpenHashMap;

public class HLA{

    public static int NEW_NODE_ADDED = 0;
    public static int HOPPING = 0;
    public static int INSERTION_NODE_ADDED = 0;
    public static int INSERTION_WITH_NO_NEW_NODE = 0;
    public static int INSERTION = 0;
    public static int READ_LENGTH = 100; //automatically gets set.
    public static double X_FACTOR = 4.0d/3.0d; //xFactor == 1 (a=4b), 4/3 (a=3b), 2 (a=2b)
    
    public static LogHandler log;
    public static boolean DEBUG;

    public static boolean OUTPUT_MERGED_MSA = false;
    public static boolean

    public HLA(String[] hlaList, String nomGFile){
	this.hlaName2Graph = new HashMap<String, HLAGraph>();
	this.hlaName2typingSequences = new HashMap<String, ArrayList<HLASequence>>();
	this.loadGraphs(hlaList, nomGFile);
    }
    

    //loads HLAGraphs as well as nomG typing sequences
    private void loadGraphs(String[] hlaList, String nomGFile){
	String tmpDir = "../db";
	HLA.log.appendln("Merging HLA sequences and building HLA graphs");
	int i;
	NomG nomG = new NomG();
	nomG.loadHlaGene2Groups(nomGFile);
	for(i=0; i<hlaList.length; i++){
	    HLA.log.appendln("processing HLA gene:\t" + hlaList[i]);
	    MergeMSFs mm = new MergeMSFs();
	    if(!mm.merge(tmpDir + File.separator +  hlaList[i] + "_nuc.txt", tmpDir + File.separator + hlaList[i] + "_gen.txt", HLA.OUTPUT_MERGED_MSA)){
		System.err.println("ERROR in MSA merging. CANNOT proceed further. Exiting..");
		System.exit(-1);
	    }
	    this.hlaName2Graph.put(hlaList[i], new HLAGraph(mm.getListOfSequences()));
	    this.hlaName2typingSequences.put(hlaList[i], mm.formDataBase(nomG.getGroups(hlaList[i])));
	    this.hlaName2Graph.get(hlaList[i]).setTypingSequences(this.hlaName2typingSequences.get(hlaList[i]));
	    this.outputTypingSequences(hlaList[i]);
	}
	HLA.log.appendln("Done building\t" + i + "\tgraphs.");
    }
    
    public void outputTypingSequences(String hgn){
	ArrayList<HLASequence> typingSeqs = this.hlaName2typingSequences.get(hgn);

	BufferedWriter bw = null;
	try{
	    bw = new BufferedWriter(new FileWriter(hgn + "_typingDB.fa"));
	    for(HLASequence h : typingSeqs){
		bw.write(Bubble.stripPadding(h.toString()));
		//bw.write(">" + h.getGroup().getGroupString() + "\n");
		//bw.write();
	    }
	    bw.close();
	}catch(IOException ioe){
	    ioe.printStackTrace();
	}
    }

    public void loadReads(File[] bams) throws IOException{
	
	int count = 0;
	int numOp = 0;
	
	for(File bam : bams){
	    HLA.log.appendln("Loading reads from:\t" + bam.getName());
	    Object2IntOpenHashMap<String> readLoadingSet = new Object2IntOpenHashMap<String>();
	    //if(pairedend){
	    //readLoadingSet = new Object2IntOpenHashMap<String>();
	    readLoadingSet.defaultReturnValue(0);
	    //}
	    
	    final SamReader reader = SamReaderFactory.makeDefault().open(bam);
	    for(final SAMRecord samRecord : reader){
		if(count == 0){
		    HLA.READ_LENGTH = samRecord.getReadLength();
		    HLA.log.appendln("Setting HLA.READ_LEGNTH = " + HLA.READ_LENGTH);
		}
		    //HLA.READ_LENGTH = (samRecord.getReadLength() > 300) ? 100 : samRecord.getReadLength();
		//System.out.println(samRecord.getCigarString());
		//samRecord
		
		//added checking to process reads matching to HLA-type sequences
		//discarding decoy hits (DQB2, DQA2)
		boolean qc = false;
		if( (samRecord.getReferenceName().indexOf("*") > -1) 
		    && !samRecord.getReadUnmappedFlag() 
		    && !samRecord.isSecondaryOrSupplementary() 
		    && !this.startWIns(samRecord)){
		    //		    && (qc==this.qcCheck(samRecord)) ){
		    count++;
		    if(samRecord.getReadPairedFlag())
			numOp += processRecord(samRecord, readLoadingSet);
		    else
			numOp += processRecordUnpaired(samRecord);
		}/*else{
		    if(!qc){
			HLA.log.appendln("SKIPPING LOW-QUAL READS");
			HLA.log.appendln(samRecord.getSAMString());
		    }
		    }*/
		if(count%10000 == 0)
		    HLA.log.appendln("Processed 10000 reads...");
	    }
	    reader.close();
	}
	HLA.log.appendln("Loaded a total of " + count + " mapped reads.");
	HLA.log.appendln("A total of " + numOp + " bases");
    }
    
    public void updateErrorProb(){
	HLA.log.appendln("------------ UPDATING error probabilities of each edge ---------");
	Iterator itr = this.hlaName2Graph.keySet().iterator();
	while(itr.hasNext()){
	    this.hlaName2Graph.get(itr.next()).updateEdgeWeightProb();
	}
	HLA.log.appendln("------------     DONE UPDATING error probabilities     ---------");
    }
    
    //assume interleaved SAMRecord
    public int processRecord(SAMRecord sr, Object2IntOpenHashMap<String> readLoadingSet){
	int totalOp = 0;
	String hlagene = HLA.extractHLAGeneName(sr.getReferenceName());
	HLAGraph hg = this.hlaName2Graph.get(hlagene);
	//hg.traverse();
	if(hg != null){
	    if(hg.isClassI()){
		boolean qc = this.qcCheck(sr);
		if(!qc)
		    return 0;
	    }
	    int readnum = readLoadingSet.getInt(sr.getReadName());
	    //no such read has been read. return value of 0 means the hashSet doesn't have the read
	    if(readnum == 0){
		readnum = sr.getFirstOfPairFlag() ? HLA.readNum : 0-HLA.readNum;
		//(132)  (-5695)  (-5137)  (5137)  (-2567)  (-1363)  (5143)
		
		readLoadingSet.put(sr.getReadName(), HLA.readNum);
		HLA.readNum++;
	    }else
		readnum = sr.getFirstOfPairFlag() ? readnum : 0-readnum;
	    
	    /*if(readnum == 2197 || readnum == -2197 || readnum == 2199 || readnum == -2199 || readnum == 2196 || readnum == -2196 || readnum ==2198 || readnum ==-2198 || readnum ==267 || readnum == -267){
		HLA.log.appendln("readnumPROBLEM( " + readnum + "):" + sr.getReadName());
		HLA.log.appendln(sr.getSAMString());
		}*/
	    /*
	    if(readnum == 1365 || readnum == -3991 || readnum == 5164 || readnum == -415 || readnum == 780 || readnum == -631 ){//if(readnum == 5137 || readnum == -5137 || readnum == 5143 || readnum == 132 || readnum == -5695 || readnum == -2567 || readnum == -1363){
	    HLA.log.appendln("readnumPROBLEM( " + readnum + "):" + sr.getReadName());
		}*/

	    totalOp += hg.addWeight(sr, readnum);//HLA.readNum);
	    //HLA.readNum++;
	}else{
	    ;//HLA.log.appendln("UNKNOWN HLA GENE: " + hlagene);
	}
	return totalOp;
    }
    
    public boolean startWIns(SAMRecord sr){
	Cigar cigar = sr.getCigar();
	if(cigar == null){
	    return true;
	}else{
	    CigarOperator op = cigar.getCigarElements().get(0).getOperator();
	    if(op == CigarOperator.I){
		if(HLA.DEBUG)
		    HLA.log.appendln("SKIPPING(Start with Insertion):\t" + sr.getReadName());
		return true;
		
	    }
	}
	return false;
    }

    public boolean qcCheck(SAMRecord sr){
	boolean debug = HLA.DEBUG;
	Cigar cigar = sr.getCigar();
	int rLen = sr.getReadLength();
	int effectiveLen = 0;
	if(cigar==null) 
	    return false;
	else{
	    for(final CigarElement ce : cigar.getCigarElements()){
		CigarOperator op = ce.getOperator();
		int cigarLen = ce.getLength();
		switch(op)
		    {
		    case M:
			{
			    effectiveLen += cigarLen;
			    break;
			}
		    case I:
			{
			    effectiveLen += cigarLen;
			    break;
			}
		    default: 
			break;
		    }
	    }
	}
	if(debug){
	    HLA.log.appendln(sr.getSAMString());
	    HLA.log.appendln("EffectiveLen:\t" + effectiveLen);
	    HLA.log.appendln("ReadLen:\t" + rLen);
	}
	Integer i = sr.getIntegerAttribute("NM");
	int nm = 0;
	if(i!=null)
	    nm = i.intValue();
	if(debug)
	    HLA.log.appendln("NM=\t" + nm);
	//if((effectiveLen*1.0d)/(rLen*1.0d) > 0.7d && (nm < 3 || ((nm*1.0d)/(effectiveLen*1.0d)) < 0.05d)){
	if(nm < 16){
	    if(debug)
		HLA.log.appendln("PASSWED QC");
	    return true;
	}
	if(debug)
	    HLA.log.appendln("FAILED QC");
	HLA.log.appendln(sr.getSAMString());
	return false;
    }


    public int processRecordUnpaired(SAMRecord sr){
	int totalOp = 0;
	String hlagene = HLA.extractHLAGeneName(sr.getReferenceName());
	HLAGraph hg = this.hlaName2Graph.get(hlagene);
	//hg.traverse();
	if(hg != null){
	    if(hg.isClassI()){
		boolean qc = this.qcCheck(sr);
		if(!qc)
		    return 0;
	    }
	    totalOp += hg.addWeight(sr, HLA.readNum);
	    HLA.readNum++;
	}
	return totalOp;
    }

    /*    
    public int processRecord(SAMRecord sr){
	int totalOp = 0;
	String hlagene = HLA.extractHLAGeneName(sr.getReferenceName());
	HLAGraph hg = this.hlaName2Graph.get(hlagene);
	//hg.traverse();
	if(hg != null){
	    totalOp += hg.addWeight(sr, HLA.readNum);
	    HLA.readNum++;
	}else{
	    ;
	}
	return totalOp;
    }
    */
    public void setNames(){
	this.hlaName2Graph.get("A").setHLAGeneName("A");
	this.hlaName2Graph.get("B").setHLAGeneName("B");
	this.hlaName2Graph.get("C").setHLAGeneName("C");
	this.hlaName2Graph.get("DQA1").setHLAGeneName("DQA1");
	this.hlaName2Graph.get("DQB1").setHLAGeneName("DQB1");
	this.hlaName2Graph.get("DRB1").setHLAGeneName("DRB1");
    }

    public void printWeights(){
	this.hlaName2Graph.get("A").traverseAndWeights();
	this.hlaName2Graph.get("B").traverseAndWeights();
	this.hlaName2Graph.get("C").traverseAndWeights();
	this.hlaName2Graph.get("DQA1").traverseAndWeights();
	this.hlaName2Graph.get("DQB1").traverseAndWeights();
	this.hlaName2Graph.get("DRB1").traverseAndWeights();
    }
    
    public void printBoundaries(){
	this.hlaName2Graph.get("A").getRefAllele().printBoundaries();
	this.hlaName2Graph.get("B").getRefAllele().printBoundaries();
	this.hlaName2Graph.get("C").getRefAllele().printBoundaries();
	this.hlaName2Graph.get("DQA1").getRefAllele().printBoundaries();
	this.hlaName2Graph.get("DQB1").getRefAllele().printBoundaries();
	this.hlaName2Graph.get("DRB1").getRefAllele().printBoundaries();
    }

    public void removeUnused(){
	this.hlaName2Graph.get("A").removeUnused();
	this.hlaName2Graph.get("B").removeUnused();
	this.hlaName2Graph.get("C").removeUnused();
	this.hlaName2Graph.get("DQA1").removeUnused();
	this.hlaName2Graph.get("DQB1").removeUnused();
	this.hlaName2Graph.get("DRB1").removeUnused();
    }

    public void flattenInsertionNodes(){
	this.hlaName2Graph.get("A").flattenInsertionNodes();
	this.hlaName2Graph.get("B").flattenInsertionNodes();
	this.hlaName2Graph.get("C").flattenInsertionNodes();
	this.hlaName2Graph.get("DQA1").flattenInsertionNodes();
	this.hlaName2Graph.get("DQB1").flattenInsertionNodes();
	this.hlaName2Graph.get("DRB1").flattenInsertionNodes();
    }


    public void printStartEndNodes(){
	this.hlaName2Graph.get("A").printStartEndNodeInfo();
	this.hlaName2Graph.get("B").printStartEndNodeInfo();
	this.hlaName2Graph.get("C").printStartEndNodeInfo();
	this.hlaName2Graph.get("DQA1").printStartEndNodeInfo();
	this.hlaName2Graph.get("DQB1").printStartEndNodeInfo();
	this.hlaName2Graph.get("DRB1").printStartEndNodeInfo();
    }

    public void countBubbles(){
	this.hlaName2Graph.get("A").countBubbles();
	this.hlaName2Graph.get("B").countBubbles();
	this.hlaName2Graph.get("C").countBubbles();
	this.hlaName2Graph.get("DQA1").countBubbles();
	this.hlaName2Graph.get("DQB1").countBubbles();
	this.hlaName2Graph.get("DRB1").countBubbles();
    }

    public void countBubblesAndMerge(StringBuffer rb){
	this.hlaName2Graph.get("A").countBubblesAndMerge(rb);
	this.hlaName2Graph.get("B").countBubblesAndMerge(rb);
	this.hlaName2Graph.get("C").countBubblesAndMerge(rb);
	this.hlaName2Graph.get("DQA1").countBubblesAndMerge(rb);
	this.hlaName2Graph.get("DQB1").countBubblesAndMerge(rb);
	this.hlaName2Graph.get("DRB1").countBubblesAndMerge(rb);
    }
    /*
    public void countBubblesAndMerge(){
	this.hlaName2Graph.get("DQA1").countBubblesAndMerge();
	this.hlaName2Graph.get("DQB1").countBubblesAndMerge();
	this.hlaName2Graph.get("DRB1").countBubblesAndMerge();
	this.hlaName2Graph.get("A").countBubblesAndMerge();
	this.hlaName2Graph.get("B").countBubblesAndMerge();
	this.hlaName2Graph.get("C").countBubblesAndMerge();
	}*/
    
    public void countStems(){
	this.hlaName2Graph.get("A").countStems();
	this.hlaName2Graph.get("B").countStems();
	this.hlaName2Graph.get("C").countStems();
	this.hlaName2Graph.get("DQA1").countStems();
	this.hlaName2Graph.get("DQB1").countStems();
	this.hlaName2Graph.get("DRB1").countStems();
    }

    public void removeStems(){
	this.hlaName2Graph.get("A").removeStems();
	this.hlaName2Graph.get("B").removeStems();
	this.hlaName2Graph.get("C").removeStems();
	this.hlaName2Graph.get("DQA1").removeStems();
	this.hlaName2Graph.get("DQB1").removeStems();
	this.hlaName2Graph.get("DRB1").removeStems();
    }

    public void setFileName(String f){
	this.outfilename = f;
	this.hlaName2Graph.get("A").setFileName(f);
	this.hlaName2Graph.get("B").setFileName(f);
	this.hlaName2Graph.get("C").setFileName(f);
	this.hlaName2Graph.get("DQA1").setFileName(f);
	this.hlaName2Graph.get("DQB1").setFileName(f);
	this.hlaName2Graph.get("DRB1").setFileName(f);
    }


    public void writeResults(StringBuffer rb){
	BufferedWriter bw = null;
	try{
	    bw = new BufferedWriter(new FileWriter(this.outfilename + ".result"));
	    bw.write(rb.toString());
	    bw.close();
	}catch(IOException ioe){
	    ioe.printStackTrace();
	}
    }
    
    public static void main(String[] args) throws IOException{
	
	if(args.length < 2){
	    System.err.println("USAGE: java -jar <PATH-TO>/Kourami.jar <bamfile1> <bamfile2> ... <bamfileN> <outfilename>");
	    System.exit(1);
	}

	String[] list = {"A" , "B" , "C" , "DQA1" , "DQB1" , "DRB1"};
	File[] bamfiles = null;
	String outfilename = null;
	if(args.length >1){
	    bamfiles = new File[args.length-1];
	    for(int i=0; i<args.length-1; i++)
		bamfiles[i]=new File(args[i]);
	    outfilename = args[args.length-1];
	}
	
	HLA.log = new LogHandler(outfilename);
	for(int i =0; i<args.length;i++)
	    HLA.log.append(args[i] + "|\t");
	HLA.log.appendln();

	try{
	    HLA hla = new HLA(list, "../db/hla_nom_g.txt");
	    //sets HLA geneNames to each graph.
	    hla.setNames();
	    //hla.printBoundaries();
	    //1. bubble counting before loading reads.
	    //System.err.println("----------------BUBBLE COUNTING: REF GRAPH--------------");
	    //HLA.log.appendln("----------------BUBBLE COUNTING: REF GRAPH--------------");
	    
	    //hla.countStems();
	    
	    System.err.println("---------------- READ LOADING --------------");
	    HLA.log.appendln("---------------- READ LOADING --------------");
	    
	    hla.loadReads(bamfiles); 
	    hla.setFileName(outfilename);
	    
	    System.err.println("---------------- GRAPH CLEANING --------------");
	    HLA.log.appendln("---------------- GRAPH CLEANING --------------");
	    	    	    
	    hla.flattenInsertionNodes();
	    hla.removeUnused();
	    hla.removeStems();
	    //hla.countStems();
	    
	    /*updating error prob*/
	    hla.updateErrorProb();
	    
	    hla.log.flush();
	    
	    StringBuffer resultBuffer = new StringBuffer();
	    
	    hla.countBubblesAndMerge(resultBuffer);
	    
	    hla.writeResults(resultBuffer);
	}catch(Exception e){
	    e.printStackTrace();
	    HLA.log.outToFile();
	    System.exit(-1);
	}
	/*printingWeights*/
	//hla.printWeights();
	HLA.log.outToFile();
	//public static int NEW_NODE_ADDED = 0;
	//public static int HOPPING = 0;
	//public static int INSERTION_NODE_ADDED = 0;
    	HLA.log.appendln("NEW_NODE_ADDED:\t" + HLA.NEW_NODE_ADDED);
	HLA.log.appendln("HOPPPING:\t" + HLA.HOPPING);
	HLA.log.appendln("INSERTION_NODE_ADDED:\t" + HLA.INSERTION_NODE_ADDED);
	HLA.log.appendln("INSERTION_WITH_NO_NEW_NODE:\t" + HLA.INSERTION_WITH_NO_NEW_NODE);
	HLA.log.appendln("INSERTION_COUNTS:\t" + HLA.INSERTION);
    }

    private static String extractHLAGeneName(String g){
	//if(g.indexOf("*") < 0)
	//return null;
	return g.substring(0,g.indexOf("*"));
    }
    
    
    public static int readNum = 1;
    private HashMap<String, HLAGraph> hlaName2Graph;
    private HashMap<String, ArrayList<HLASequence>> hlaName2typingSequences;
    private String outfilename;
}
