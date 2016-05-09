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
	    totalOp += hg.addWeight(sr);
	}else{
	    ;//System.err.println("UNKNOWN HLA GENE: " + hlagene);
	}
	return totalOp;
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

    public static void main(String[] args) throws IOException{
	String[] list = {"A" , "B" , "C" , "DQA1" , "DQB1" , "DRB1"};
	//list[0] = args[1];
	if(args.length > 1){
	    list = new String[1];
	    list[0] = args[1];
	}
	HLA hla = new HLA(list);
	//hla.printBoundaries();
	hla.loadReads(new File(args[0]));
	hla.updateErrorProb();
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
    
    

    private HashMap<String, HLAGraph> hlaName2Graph;
}
