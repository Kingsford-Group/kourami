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
	for(int i=0; i<hlaList.length; i++){
	
	}
    }
    
    public void loadReads(File bam) throws IOException{
	final SamReader reader = SamReaderFactory.makeDefault().open(bam);
	for(final SAMRecord samRecord : reader){
	    //System.out.println(samRecord.getCigarString());
	    //samRecord
	    if(samRecord.getReadUnmappedFlag())
		processRecord(samRecord);
	}
	reader.close();
    }
    
    public void processRecord(SAMRecord sr){
	String hlagene = HLA.extractHLAGeneName(sr.getReferenceName());
	HLAGraph hg = this.hlaName2Graph.get(hlagene);
	hg.addWeight(sr);
    }

    public static void main(String[] args) throws IOException{
	new HLA().loadReads(new File(args[0]));
    }

    private static String extractHLAGeneName(String g){
	return g.substring(0,g.indexOf("*"));
    }
    
    

    private HashMap<String, HLAGraph> hlaName2Graph;
}
