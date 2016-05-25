import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

import java.io.File;
import java.io.IOException;
import java.util.HashMap;
import java.util.Iterator;

public class HLA{

    public static int NEW_NODE_ADDED = 0;
    public static int HOPPING = 0;
    public static int INSERTION_NODE_ADDED = 0;
    public static int INSERTION_WITH_NO_NEW_NODE = 0;
    public static int INSERTION = 0;
    
    public HLA(String[] hlaList){
	this.hlaName2Graph = new HashMap<String, HLAGraph>();
	this.loadGraphs(hlaList);
    }
    
    private void loadGraphs(String[] hlaList){
	String tmpDir = "/home/heewookl/utilities/msfs/";
	System.err.println("Merging HLA sequences and building HLA graphs");
	int i;
	for(i=0; i<hlaList.length; i++){
	    System.err.println("processing HLA gene:\t" + hlaList[i]);
	    MergeMSFs mm = new MergeMSFs();
	    mm.merge(tmpDir + hlaList[i] + "_nuc.txt", tmpDir + hlaList[i] + "_gen.txt");
	    //mm.merge(hlaList[i] + "_nuc_merged.txt", hlaList[i] + "_gen_merged.txt");
	    //mm.merge(hlaList[i] + "_nuc_short_test.txt", hlaList[i] + "_gen_short_test.txt");
	    //mm.merge(hlaList[i] + "_nuc_long_test.txt", hlaList[i] + "_gen_long_test.txt");
	    this.hlaName2Graph.put(hlaList[i], new HLAGraph(mm.getListOfSequences()));
	}
	System.err.println("Done building\t" + i + "\tgraphs.");
    }
    
    public void loadReads(File bam) throws IOException{
	System.err.println("Loading reads from:\t" + bam.getName());
	int count = 0;
	int numOp = 0;
	final SamReader reader = SamReaderFactory.makeDefault().open(bam);
	for(final SAMRecord samRecord : reader){
	    //System.out.println(samRecord.getCigarString());
	    //samRecord
	    if(!samRecord.getReadUnmappedFlag()){
		count++;
		numOp += processRecord(samRecord);
	    }
	    if(count%10000 == 0)
		System.err.println("Processed 10000 reads...");
	}
	reader.close();
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
    
    public int processRecord(SAMRecord sr){
	int totalOp = 0;
	String hlagene = HLA.extractHLAGeneName(sr.getReferenceName());
	HLAGraph hg = this.hlaName2Graph.get(hlagene);
	//hg.traverse();
	if(hg != null){
	    totalOp += hg.addWeight(sr, HLA.readNum);
	    HLA.readNum++;
	}else{
	    ;//System.err.println("UNKNOWN HLA GENE: " + hlagene);
	}
	return totalOp;
    }

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

    public static void main(String[] args) throws IOException{
	String[] list = {"A" , "B" , "C" , "DQA1" , "DQB1" , "DRB1"};
	//list[0] = args[1];
	if(args.length > 1){
	    list = new String[1];
	    list[0] = args[1];
	}
	HLA hla = new HLA(list);
	//sets HLA geneNames to each graph.
	hla.setNames();
	//hla.printBoundaries();
	//1. bubble counting before loading reads.
	System.err.println("----------------BUBBLE COUNTING: REF GRAPH--------------");
	
	hla.countStems();
	hla.loadReads(new File(args[0]));

	//2. bubble counting after loading reads
	System.err.println("----------------BUBBLE COUNTING: POST-read loading--------------");
	//hla.countBubbles();
	//hla.countStems();
	//hla.removeStems();
	//	hla.countStems();

	hla.printBoundaries();

	hla.printStartEndNodes();

	//hla.countBubbles();
	
	hla.removeUnused();
	
	
	//hla.countStems();
	
	//hla.removeStems();
	//hla.countStems();
	
	hla.flattenInsertionNodes();
	
	hla.removeStems();
	hla.countStems();
	
	/*updating error prob*/
	hla.updateErrorProb();
	
	hla.countBubbles();

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
	return g.substring(0,g.indexOf("*"));
    }
    
    
    public static int readNum = 0;
    private HashMap<String, HLAGraph> hlaName2Graph;
}
