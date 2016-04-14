import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

import java.io.File;
import java.io.IOException;
import java.util.HashMap;

public class HLA{

    public HLA(String[] hlaList){
	this.hlaName2Graph = new HashMap<String, HLAGraph>();
	this.loadGraphs(hlaList);
    }
    
    private void loadGraphs(String[] hlaList){
	System.err.println("Merging HLA sequences and building HLA graphs");
	int i;
	for(i=0; i<hlaList.length; i++){
	    System.err.println("processing HLA gene:\t" + hlaList[i]);
	    MergeMSFs mm = new MergeMSFs();
	    mm.merge(hlaList[i] + "_nuc_merged.txt", hlaList[i] + "_gen_merged.txt");
	    //mm.merge(hlaList[i] + "_nuc_short_test.txt", hlaList[i] + "_gen_short_test.txt");
	    //mm.merge(hlaList[i] + "_nuc_long_test.txt", hlaList[i] + "_gen_long_test.txt");
	    this.hlaName2Graph.put(hlaList[i], new HLAGraph(mm.getListOfSequences()));
	}
	System.err.println("Done building\t" + i + "\tgraphs.");
    }
    
    public void loadReads(File bam) throws IOException{
	System.err.println("Loading reads from:\t" + bam.getName());
	int count = 0;
	final SamReader reader = SamReaderFactory.makeDefault().open(bam);
	for(final SAMRecord samRecord : reader){
	    //System.out.println(samRecord.getCigarString());
	    //samRecord
	    if(!samRecord.getReadUnmappedFlag()){
		count++;
		processRecord(samRecord);
	    }
	    if(count%20 == 0)
		System.err.println("Processed 20 reads.");
	}
	reader.close();
	System.err.println("Loaded a total of " + count + " mapped reads.");
    }
    
    public void processRecord(SAMRecord sr){
	String hlagene = HLA.extractHLAGeneName(sr.getReferenceName());
	HLAGraph hg = this.hlaName2Graph.get(hlagene);
	//hg.traverse();
	if(hg != null){
	    hg.addWeight(sr);
	}else{
	    System.err.println("UNKNOWN HLA GENE: " + hlagene);
	}
    }

    public static void main(String[] args) throws IOException{
	String[] list = new String[1];
	list[0] = args[1];
	new HLA(list).loadReads(new File(args[0]));
    }

    private static String extractHLAGeneName(String g){
	return g.substring(0,g.indexOf("*"));
    }
    
    

    private HashMap<String, HLAGraph> hlaName2Graph;
}
