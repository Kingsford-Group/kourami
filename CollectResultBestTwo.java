import java.io.*;
import java.util.*;

/*
 * This version uses BestTwo.
 *
 *
 */

public class CollectResultBestTwo{
    
    private HashMap<String, String> hlaName2Typing;
    
    //private HashMap<String, SortedRecords> hlaLocus2Alleles;
    private HashMap<String, BestTwo> hlaLocus2Alleles;
    
    public static void main(String[] args){
	CollectResultBestTwo obj = new CollectResultBestTwo();
	obj.loadNames(args[0]); //names_file
	obj.processQuant(args[1]); //quant file
	obj.output();
	//obj.outToFile(args[2]);
    }
    
    //public SortedRecords getBestTwoFromHLAType(String typeStr){
    public BestTwo getBestTwoFromHLAType(String typeStr){
	//System.out.println("TypeStr\t" + typeStr + "\tKey is:" + typeStr.substring(0, typeStr.indexOf("*")));
	
	return this.hlaLocus2Alleles.get(typeStr.substring(0, typeStr.indexOf("*")));
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


    /*    public void outToFile(String f){
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
	}*/
    
    
    public void processQuant(String qf){
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
		    //SortedRecords curb2 = this.getBestTwoFromHLAType(this.hlaName2Typing.get(tokens[0]));
		    BestTwo curb2 = this.getBestTwoFromHLAType(this.hlaName2Typing.get(tokens[0]));
		    curb2.update(this.hlaName2Typing.get(tokens[0]), tpm);
		}
	    }
	    br.close();
	}catch(IOException ioe){
	    ioe.printStackTrace();
	}
    }
    
    public void loadNames(String namesF){
	BufferedReader br = null;
	try{
	    br = new BufferedReader(new FileReader(namesF));
	    String curline = "";
	    while((curline=br.readLine())!=null){
		String[] tokens = curline.substring(1).split("\\s+");
		String hlaname = tokens[0];
		String type = tokens[1];
		//System.out.println("Putting [" + hlaname+"],["+ type +"]");
		this.hlaName2Typing.put(hlaname, type);
	    }
	    br.close();
	}catch(IOException ioe){
	    ioe.printStackTrace();
	}
    }
    
    public CollectResultBestTwo(){
	this.hlaName2Typing = new HashMap<String, String>();
	this.hlaLocus2Alleles = new HashMap<String, BestTwo>();
	this.hlaLocus2Alleles.put("A", new BestTwo());
	this.hlaLocus2Alleles.put("B", new BestTwo());
	this.hlaLocus2Alleles.put("C", new BestTwo());
	this.hlaLocus2Alleles.put("DQA1", new BestTwo());
	this.hlaLocus2Alleles.put("DMA", new BestTwo());
	this.hlaLocus2Alleles.put("DMB", new BestTwo());
	this.hlaLocus2Alleles.put("DOA", new BestTwo());
	this.hlaLocus2Alleles.put("DOB", new BestTwo());
	this.hlaLocus2Alleles.put("DPA1", new BestTwo());
	this.hlaLocus2Alleles.put("DPB1", new BestTwo());
	this.hlaLocus2Alleles.put("DQA1", new BestTwo());
	this.hlaLocus2Alleles.put("DQB1", new BestTwo());
	this.hlaLocus2Alleles.put("DRA", new BestTwo());
	this.hlaLocus2Alleles.put("DRB1", new BestTwo());
	this.hlaLocus2Alleles.put("DRB3", new BestTwo());
	this.hlaLocus2Alleles.put("DRB4", new BestTwo());
	this.hlaLocus2Alleles.put("E", new BestTwo());
	this.hlaLocus2Alleles.put("F", new BestTwo());
	this.hlaLocus2Alleles.put("G", new BestTwo());
	this.hlaLocus2Alleles.put("H", new BestTwo());
	this.hlaLocus2Alleles.put("J", new BestTwo());
	this.hlaLocus2Alleles.put("K", new BestTwo());
	this.hlaLocus2Alleles.put("L", new BestTwo());
	this.hlaLocus2Alleles.put("MICA", new BestTwo());
	this.hlaLocus2Alleles.put("MICB", new BestTwo());
	this.hlaLocus2Alleles.put("P", new BestTwo());
	this.hlaLocus2Alleles.put("TAP1", new BestTwo());
	this.hlaLocus2Alleles.put("TAP2", new BestTwo());
	this.hlaLocus2Alleles.put("V", new BestTwo());
    }

}


class BestTwo{
    
    double tpm1;
    double tpm2;
    StringBuffer name1; // types
    StringBuffer name2;

    public String extractUptoXdigit(int n){
	//System.out.print("-->\t" + name1 + "\t" + name2);
	int len1 = this.name1.length();
	int len2 = this.name2.length();
	StringBuffer out = new StringBuffer("-->\t");
	if(len1 > 0){
	    int starPos1 = name1.indexOf("*")+1;
	    out.append(name1.substring(0,name1.indexOf("*")) + "*" + name1.substring(starPos1, starPos1+n+((n-1)/2)));
	    if(len2 > 0){
		int starPos2 = name2.indexOf("*")+1;
		out.append("\t" + name2.substring(0,name2.indexOf("*")) + "*" + name2.substring(starPos2, starPos2+n+((n-1)/2)));
	    }
	}else
	    out.append("NO ALLELE DETECTED");
	return out.toString();
    }
    
    public BestTwo(){
	this.tpm1 = 0.0d;
	this.tpm2 = 0.0d;
	this.name1 = new StringBuffer();
	this.name2 = new StringBuffer();
    }
    
    public void update(String name, double tpm){
	if(tpm > this.tpm2){
	    if(tpm <= this.tpm1){
		this.tpm2 = tpm;
		this.name2 = new StringBuffer(name);
	    }else{
		this.tpm2 = this.tpm1;
		this.name2 = this.name1;
		this.tpm1 = tpm;
		this.name1 = new StringBuffer(name);
	    }
	}
	
    }
}
