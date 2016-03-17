import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

import java.io.File;
import java.io.IOException;
import java.util.HashMap;

public class HLA{

    public HLA(){
	this.loadGraphs();
    }
    
    private void loadGraphs(){
	this.hlaName2Graph = null;
    }
    
    public void loadReads(File bam) throws IOException{
	final SamReader reader = SamReaderFactory.makeDefault().open(bam);
	for(final SAMRecord samRecord : reader){
	    //System.out.println(samRecord.getCigarString());
	    //samRecord
	    if(samRecord.getReadUnmappedFlag())
		process(samRecord);
	}
	reader.close();
    }
    
    public void processRecord(SAMRecord sr){
	String hlagene = HLA.extractHLAGeneName(sr.getReferenceName());
	HLAGraph g = this.hlaName2Graph.get(hlagene);
	
    }

    public static void main(String[] args) throws IOException{
	new HLA().loadReads(new File(args[0]));
    }

    private static String extractHLAGeneName(String g){
	return g.substring(0,g.indexOf("*"));
    }
    
    

    private HashMap<String, HLAGraph> hlaName2Graph;
}
