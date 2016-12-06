import htsjdk.samtools.SAMRecord;
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
    public static int READ_LENGTH = 100;
    public static double X_FACTOR = 4.0d/3.0d; //xFactor == 1 (a=4b), 4/3 (a=3b), 2 (a=2b)
    
    public HLA(String[] hlaList, String nomGFile){
	this.hlaName2Graph = new HashMap<String, HLAGraph>();
	this.hlaName2typingSequences = new HashMap<String, ArrayList<HLASequence>>();
	this.loadGraphs(hlaList, nomGFile);
    }
    

    //loads HLAGraphs as well as nomG typing sequences
    private void loadGraphs(String[] hlaList, String nomGFile){
	String tmpDir = "/home/heewookl/utilities/msfs/";
	//String tmpDir = "/home/heewookl/utilities/msfs/WithAnswersOutNA12878/";
	System.err.println("Merging HLA sequences and building HLA graphs");
	int i;
	NomG nomG = new NomG();
	nomG.loadHlaGene2Groups(nomGFile);
	for(i=0; i<hlaList.length; i++){
	    System.err.println("processing HLA gene:\t" + hlaList[i]);
	    MergeMSFs mm = new MergeMSFs();
	    mm.merge(tmpDir + hlaList[i] + "_nuc.txt", tmpDir + hlaList[i] + "_gen.txt");
	    //mm.merge(hlaList[i] + "_nuc_merged.txt", hlaList[i] + "_gen_merged.txt");
	    //mm.merge(hlaList[i] + "_nuc_short_test.txt", hlaList[i] + "_gen_short_test.txt");
	    //mm.merge(hlaList[i] + "_nuc_long_test.txt", hlaList[i] + "_gen_long_test.txt");
	    this.hlaName2Graph.put(hlaList[i], new HLAGraph(mm.getListOfSequences()));
	    this.hlaName2typingSequences.put(hlaList[i], mm.formDataBase(nomG.getGroups(hlaList[i])));
	    this.hlaName2Graph.get(hlaList[i]).setTypingSequences(this.hlaName2typingSequences.get(hlaList[i]));
	    this.outputTypingSequences(hlaList[i]);
	}
	System.err.println("Done building\t" + i + "\tgraphs.");
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
	    System.err.println("Loading reads from:\t" + bam.getName());
	    Object2IntOpenHashMap<String> readLoadingSet = new Object2IntOpenHashMap<String>();
	    //if(pairedend){
	    //readLoadingSet = new Object2IntOpenHashMap<String>();
	    readLoadingSet.defaultReturnValue(0);
	    //}
	    
	    final SamReader reader = SamReaderFactory.makeDefault().open(bam);
	    for(final SAMRecord samRecord : reader){
		if(count == 0){
		    HLA.READ_LENGTH = samRecord.getReadLength();
		    System.err.println("Setting HLA.READ_LEGNTH = " + HLA.READ_LENGTH);
		}
		    //HLA.READ_LENGTH = (samRecord.getReadLength() > 300) ? 100 : samRecord.getReadLength();
		//System.out.println(samRecord.getCigarString());
		//samRecord
		
		//added checking to process reads matching to HLA-type sequences
		//discarding decoy hits (DQB2, DQA2)
		if( (samRecord.getReferenceName().indexOf("*") > -1) 
		    && !samRecord.getReadUnmappedFlag() ){
		    count++;
		    if(samRecord.getReadPairedFlag())
			numOp += processRecord(samRecord, readLoadingSet);
		    else
			numOp += processRecordUnpaired(samRecord);
		}
		if(count%10000 == 0)
		    System.err.println("Processed 10000 reads...");
	    }
	    reader.close();
	}
	System.err.println("Loaded a total of " + count + " mapped reads.");
	System.err.println("A total of " + numOp + " bases");
    }
    
    public void updateErrorProb(){
	System.err.println("------------ UPDATING error probabilities of each edge ---------");
	Iterator itr = this.hlaName2Graph.keySet().iterator();
	while(itr.hasNext()){
	    this.hlaName2Graph.get(itr.next()).updateEdgeWeightProb();
	}
	System.err.println("------------     DONE UPDATING error probabilities     ---------");
    }
    
    //assume interleaved SAMRecord
    public int processRecord(SAMRecord sr, Object2IntOpenHashMap<String> readLoadingSet){
	int totalOp = 0;
	String hlagene = HLA.extractHLAGeneName(sr.getReferenceName());
	HLAGraph hg = this.hlaName2Graph.get(hlagene);
	//hg.traverse();
	if(hg != null){
	    int readnum = readLoadingSet.getInt(sr.getReadName());
	    //no such read has been read. return value of 0 means the hashSet doesn't have the read
	    if(readnum == 0){
		readnum = sr.getFirstOfPairFlag() ? HLA.readNum : 0-HLA.readNum;
		//(132)  (-5695)  (-5137)  (5137)  (-2567)  (-1363)  (5143)
		
		readLoadingSet.put(sr.getReadName(), HLA.readNum);
		HLA.readNum++;
	    }else
		readnum = sr.getFirstOfPairFlag() ? readnum : 0-readnum;
	    
	    if(readnum == 1365 || readnum == -3991 || readnum == 5164 || readnum == -415 || readnum == 780 || readnum == -631 ){//if(readnum == 5137 || readnum == -5137 || readnum == 5143 || readnum == 132 || readnum == -5695 || readnum == -2567 || readnum == -1363){
		System.err.println("readnumPROBLEM( " + readnum + "):" + sr.getReadName());
	    }

	    totalOp += hg.addWeight(sr, readnum);//HLA.readNum);
	    //HLA.readNum++;
	}else{
	    ;//System.err.println("UNKNOWN HLA GENE: " + hlagene);
	}
	return totalOp;
    }

    public int processRecordUnpaired(SAMRecord sr){
	int totalOp = 0;
	String hlagene = HLA.extractHLAGeneName(sr.getReferenceName());
	HLAGraph hg = this.hlaName2Graph.get(hlagene);
	//hg.traverse();
	if(hg != null){
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
	for(int i =0; i<args.length;i++)
	    System.err.print(args[i] + "|\t");
	System.err.println();
	String[] list = {"A" , "B" , "C" , "DQA1" , "DQB1" , "DRB1"};
	//list[0] = args[1];
	//boolean pairedend = true;
	//if(args.length > 2){
	//list = new String[1];
	//    list[0] = args[2];
	    //pairedend = ( (args[2].equals("N") || args[2].equals("n")) ? false : true);
	//}
	File[] bamfiles = null;
	String outfilename = null;
	if(args.length >1){
	    bamfiles = new File[args.length-1];
	    for(int i=0; i<args.length-1; i++)
		bamfiles[i]=new File(args[i]);
	    outfilename = args[args.length-1];
	}
	
	HLA hla = new HLA(list, "/home/heewookl/utilities/hla_nom_g.txt");
	//HLA hla = new HLA(list, "hla_nom_g.txt");
	//HLA hla = new HLA(list, "/home/heewookl/utilities/msfs/WithAnswersOutNA12878/hla_nom_g.txt");
	//sets HLA geneNames to each graph.
	hla.setNames();
	//hla.printBoundaries();
	//1. bubble counting before loading reads.
	System.err.println("----------------BUBBLE COUNTING: REF GRAPH--------------");
	
	hla.countStems();
	hla.loadReads(bamfiles);//new File(args[0]));//, pairedend);
	
	//if(args.length > 1){
	hla.setFileName(outfilename);//args[1]);
	    //}
	

	//2. bubble counting after loading reads
	System.err.println("----------------BUBBLE COUNTING: POST-read loading--------------");
	//hla.countBubbles();
	//hla.countStems();
	//hla.removeStems();
	//	hla.countStems();

	hla.printBoundaries();

	hla.printStartEndNodes();

	//hla.countBubbles();

	hla.flattenInsertionNodes();
	
	hla.removeUnused();
	
	
	//hla.countStems();
	
	//hla.removeStems();
	//hla.countStems();
	

	
	hla.removeStems();
	hla.countStems();
	
	/*updating error prob*/
	hla.updateErrorProb();
	
	StringBuffer resultBuffer = new StringBuffer();
	hla.countBubblesAndMerge(resultBuffer);
	
	hla.writeResults(resultBuffer);

	/*printingWeights*/
	hla.printWeights();

	//public static int NEW_NODE_ADDED = 0;
	//public static int HOPPING = 0;
	//public static int INSERTION_NODE_ADDED = 0;
    	System.err.println("NEW_NODE_ADDED:\t" + HLA.NEW_NODE_ADDED);
	System.err.println("HOPPPING:\t" + HLA.HOPPING);
	System.err.println("INSERTION_NODE_ADDED:\t" + HLA.INSERTION_NODE_ADDED);
	System.err.println("INSERTION_WITH_NO_NEW_NODE:\t" + HLA.INSERTION_WITH_NO_NEW_NODE);
	System.err.println("INSERTION_COUNTS:\t" + HLA.INSERTION);
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
