import java.io.*;
import java.util.*;

public class CollectResultMulti{
    
    private HashMap<String, String> hlaName2Typing;
    
    private HashMap<String, SortedRecords> hlaLocus2Alleles;

    //private HashMap<String, Solution> hlaSolutions;//key is geneString --> 'A', 'B', ...
    
    /*public static void main(String[] args){
	CollectResultMulti obj = new CollectResultMulti();
	//obj.loadNames(args[0]); //names_file
	//obj.loadSolutions(args[1]); //solution_file
	obj.processQuant(args[1]); //quant file
	obj.output();
	obj.outToFile(args[2]);
    }
    */    
    public SortedRecords getBestTwoFromHLAType(String typeStr){
	//System.out.println("TypeStr\t" + typeStr + "\tKey is:" + typeStr.substring(0, typeStr.indexOf("*")));
	
	return this.hlaLocus2Alleles.get(typeStr.substring(0, typeStr.indexOf("*")));
    }
    
    public HashMap<String, SortedRecords> getSortedRecordsHash(){
	return this.hlaLocus2Alleles;
    }
    /*
    public double getAccuracy(String curSampleID){
	int scores[] = new int[6];
	//for each gene
	for(int i=0;i<scores.length; i++){
	    Solution s = this.hlaSolutions.get(this.index2Gene(i));
	    scores[i] = s.matchScore(hlaLocus2Alleles.get(this.index2Gene(i)), curSampleID);
	}
	}*/
    
    public String index2Gene(int i){
	if(i == 0)
	    return "A";
	else if(i == 1)
	    return "B";
	else if(i == 2)
	    return "C";
	else if(i == 3)
	    return "DQA1";
	else if(i == 4)
	    return "DQB1";
	else if(i == 5)
	    return "DRB1";
	else
	    return null;
    }
    
    public void output(){
	//Set<String> keys = this.hlaLocus2Alleles.keySet();
	//for(String s : keys){
	    
	//    System.out.println(s + "\t" + this.hlaLocus2Alleles.get(s).extractUptoXdigit(4));
	//}
	System.out.println("A   " + "\t" + this.hlaLocus2Alleles.get("A").extractUptoXdigit(4));
	System.out.println("B   " + "\t" + this.hlaLocus2Alleles.get("B").extractUptoXdigit(4));
	System.out.println("C   " + "\t" + this.hlaLocus2Alleles.get("C").extractUptoXdigit(4));
	System.out.println("DQA1" + "\t" + this.hlaLocus2Alleles.get("DQA1").extractUptoXdigit(4));
	System.out.println("DQB1" + "\t" + this.hlaLocus2Alleles.get("DQB1").extractUptoXdigit(4));
	System.out.println("DRB1" + "\t" + this.hlaLocus2Alleles.get("DRB1").extractUptoXdigit(4));
    }

    
    public void outToFile(String f){
	BufferedWriter bw = null;
	try{
	    bw = new BufferedWriter(new FileWriter(f));
	    bw.write("A   " + "\t" + this.hlaLocus2Alleles.get("A").extractUptoXdigitFull(4) + "\n");
	    bw.write("B   " + "\t" + this.hlaLocus2Alleles.get("B").extractUptoXdigitFull(4) + "\n");
	    bw.write("C   " + "\t" + this.hlaLocus2Alleles.get("C").extractUptoXdigitFull(4) + "\n");
	    bw.write("DQA1" + "\t" + this.hlaLocus2Alleles.get("DQA1").extractUptoXdigitFull(4) + "\n");
	    bw.write("DQB1" + "\t" + this.hlaLocus2Alleles.get("DQB1").extractUptoXdigitFull(4) + "\n");
	    bw.write("DRB1" + "\t" + this.hlaLocus2Alleles.get("DRB1").extractUptoXdigitFull(4) + "\n");
	    bw.close();
	}catch(IOException ioe){
	    ioe.printStackTrace();
	}
    }
    
    
    public void processQuant(String qf){
	this.initAlleles();
	BufferedReader br = null;
	try{
	    br = new BufferedReader(new FileReader(qf));
	    String curline = "";
	    while((curline = br.readLine())!=null){
		if(!curline.startsWith("Name")){
		    String[] tokens = curline.split("\\t");
		    double tpm = Double.parseDouble(tokens[3]);
		    //if(tpm > 0)
			//System.err.println(this.hlaName2Typing.get(tokens[0]) + "\t"+ curline);
		    //System.out.println("Typing : " + this.hlaName2Typing.get(tokens[0]) + "\tHLA_NAME: " + tokens[0]);
		    SortedRecords curb2 = this.getBestTwoFromHLAType(this.hlaName2Typing.get(tokens[0]));
		    curb2.update(this.hlaName2Typing.get(tokens[0]), tpm);
		}
	    }
	    br.close();
	}catch(IOException ioe){
	    ioe.printStackTrace();
	}
    }
    
    

    
    
    public void initAlleles(){
	this.hlaLocus2Alleles = new HashMap<String, SortedRecords>();
	this.hlaLocus2Alleles.put("A", new SortedRecords());
	this.hlaLocus2Alleles.put("B", new SortedRecords());
	this.hlaLocus2Alleles.put("C", new SortedRecords());
	this.hlaLocus2Alleles.put("DQA1", new SortedRecords());
	this.hlaLocus2Alleles.put("DMA", new SortedRecords());
	this.hlaLocus2Alleles.put("DMB", new SortedRecords());
	this.hlaLocus2Alleles.put("DOA", new SortedRecords());
	this.hlaLocus2Alleles.put("DOB", new SortedRecords());
	this.hlaLocus2Alleles.put("DPA1", new SortedRecords());
	this.hlaLocus2Alleles.put("DPB1", new SortedRecords());
	this.hlaLocus2Alleles.put("DQA1", new SortedRecords());
	this.hlaLocus2Alleles.put("DQB1", new SortedRecords());
	this.hlaLocus2Alleles.put("DRA", new SortedRecords());
	this.hlaLocus2Alleles.put("DRB1", new SortedRecords());
	this.hlaLocus2Alleles.put("DRB3", new SortedRecords());
	this.hlaLocus2Alleles.put("DRB4", new SortedRecords());
	this.hlaLocus2Alleles.put("E", new SortedRecords());
	this.hlaLocus2Alleles.put("F", new SortedRecords());
	this.hlaLocus2Alleles.put("G", new SortedRecords());
	this.hlaLocus2Alleles.put("H", new SortedRecords());
	this.hlaLocus2Alleles.put("J", new SortedRecords());
	this.hlaLocus2Alleles.put("K", new SortedRecords());
	this.hlaLocus2Alleles.put("L", new SortedRecords());
	this.hlaLocus2Alleles.put("MICA", new SortedRecords());
	this.hlaLocus2Alleles.put("MICB", new SortedRecords());
	this.hlaLocus2Alleles.put("P", new SortedRecords());
	this.hlaLocus2Alleles.put("TAP1", new SortedRecords());
	this.hlaLocus2Alleles.put("TAP2", new SortedRecords());
	this.hlaLocus2Alleles.put("V", new SortedRecords());
    }

    public CollectResultMulti(){
	this.hlaName2Typing = new HashMap<String, String>();
	//this.hlaSolutions = new HashMap<String, Solution>();
    }

    public CollectResultMulti(HashMap<String, String> hlaN2T, String quant){
	this.hlaName2Typing = hlaN2T;
	this.initAlleles();
	this.processQuant(quant);
	//this.hlaSolutions = new HashMap<String, Solution>();
    }
    

}

