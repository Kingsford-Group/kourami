import java.util.*;
import java.io.*;

public class Msf2MultiFasta{

    public Msf2MultiFasta(){
	hlas = new HashMap<String, Sequence>();
	keys = new ArrayList<String>();
    }
    
    public static void main(String[] args){
	if(args[0].equals("A_nuc.msf"))
	    new Msf2MultiFasta().processA(args[0]);
	if(args[0].equals("B_nuc.msf"))
	    new Msf2MultiFasta().processB(args[0]);
	if(args[0].equals("C_nuc.msf"))
	    new Msf2MultiFasta().processC(args[0]);
	if(args[0].equals("DQA1_nuc.msf"))
	    new Msf2MultiFasta().processDQA1(args[0]);
    }
    
    public void output(){
	for(int i=0; i<this.keys.size(); i++){
	    System.out.print(this.hlas.get(keys.get(i)).getFastaStr());
	}
    }

    public void processA(String f){
	BufferedReader br = null;
	try{
	    br = new BufferedReader(new FileReader(f));
	    String curline = null;
	    boolean inSequence = false;
	    int len = 0;
	    while((curline=br.readLine())!=null){
		curline = curline.trim();
		/*if(curline.contains("A*80:03")){
		    System.err.println(curline);
		    String tmp = br.readLine();//System.err.println(br.readLine());
		    if(tmp.equals("//"))
			System.err.println("equals");
		    else
			System.err.println("Not euqls:" + tmp);
			}*/
		if(inSequence){//if we are in MSA
		    //System.err.println("inSequence");
		    if(curline.equals(""))// blank line means an end of a block
			len++;
		    else{//in MSA block
			//System.err.println("len\t-->\t" + len);
			if(len > 0){//first block we don't include since we are only targetting for exon2 and exon3
			    String curName = curline.substring(0,curline.indexOf("  "));
			    String curSeq = null;
			    //System.err.println(curName);
	    
			    if(len < 15){
				//second block contains the start of exon2
				if(len == 1){
				    curSeq = curline.substring(curline.indexOf("  ")+27).replaceAll("\\s+", "");
				    //System.err.println(curline + "\n" + curSeq);
				}else if(len == 14){ // 15th block contains the end of exon3
				    curSeq = curline.substring(curline.indexOf("  ")).substring(0, 30).replaceAll("\\s+", "");
				    //System.err.println(curline + "\n" + curSeq + "|" + curSeq.length());
				}
				else
				    curSeq = curline.substring(curline.indexOf("  ")).replaceAll("\\s+","");
				if(hlas.get(curName) == null){
				    if(len != 1)
					System.err.println("Something is not right");
				    hlas.put(curName, new Sequence(curName, curSeq));
				    keys.add(curName);
				}else
				    hlas.get(curName).append(curSeq);
			    }
			}
		    }
		}else{//header section
		    //System.err.println(curline);
		    if(curline.equals("//")){
			br.readLine();
			inSequence = true;
			//System.err.println("HERE" + curline);
		    }
		}
	    }
	    br.close();
	}catch(IOException ioe){
	    ioe.printStackTrace();
	}
	this.output();
    }

    public void processB(String f){
	BufferedReader br = null;
	try{
	    br = new BufferedReader(new FileReader(f));
	    String curline = null;
	    boolean inSequence = false;
	    int len = 0;
	    while((curline=br.readLine())!=null){
		curline = curline.trim();
		/*if(curline.contains("A*80:03")){
		    System.err.println(curline);
		    String tmp = br.readLine();//System.err.println(br.readLine());
		    if(tmp.equals("//"))
			System.err.println("equals");
		    else
			System.err.println("Not euqls:" + tmp);
			}*/
		if(inSequence){//if we are in MSA
		    //System.err.println("inSequence");
		    if(curline.equals(""))// blank line means an end of a block
			len++;
		    else{//in MSA block
			//System.err.println("len\t-->\t" + len);
			if(len > 2){//first block we don't include since we are only targetting for exon2 and exon3
			    String curName = curline.substring(0,curline.indexOf("  "));
			    String curSeq = null;
			    //System.err.println(curName);
	    
			    if(len < 16){
				//fourth block contains the start of exon2
				if(len == 3){
				    curSeq = curline.substring(curline.indexOf("  ")+48).replaceAll("\\s+", "");
				    //System.err.println(curline + "\n" + curSeq);
				}else if(len == 15){ // 16th block contains the end of exon3
				    //System.err.println(15);
				    curSeq = curline.substring(curline.indexOf("  ")).substring(0, 32).replaceAll("\\s+", "");
				    //System.err.println(curline + "\n" + curSeq + "|" + curSeq.length());
				}
				else
				    curSeq = curline.substring(curline.indexOf("  ")).replaceAll("\\s+","");
				if(hlas.get(curName) == null){
				    if(len != 3)
					System.err.println("b:Something is not right");
				    hlas.put(curName, new Sequence(curName, curSeq));
				    keys.add(curName);
				}else
				    hlas.get(curName).append(curSeq);
			    }
			}
		    }
		}else{//header section
		    //System.err.println(curline);
		    if(curline.equals("//")){
			br.readLine();
			inSequence = true;
			//System.err.println("HERE" + curline);
		    }
		}
	    }
	    br.close();
	}catch(IOException ioe){
	    ioe.printStackTrace();
	}
	//this.output();
    }

    public void processC(String f){
	BufferedReader br = null;
	try{
	    br = new BufferedReader(new FileReader(f));
	    String curline = null;
	    boolean inSequence = false;
	    int len = 0;
	    while((curline=br.readLine())!=null){
		curline = curline.trim();
		/*if(curline.contains("A*80:03")){
		    System.err.println(curline);
		    String tmp = br.readLine();//System.err.println(br.readLine());
		    if(tmp.equals("//"))
			System.err.println("equals");
		    else
			System.err.println("Not euqls:" + tmp);
			}*/
		if(inSequence){//if we are in MSA
		    //System.err.println("inSequence");
		    if(curline.equals(""))// blank line means an end of a block
			len++;
		    else{//in MSA block
			//System.err.println("len\t-->\t" + len);
			if(len > 0){//first block we don't include since we are only targetting for exon2 and exon3
			    String curName = curline.substring(0,curline.indexOf("  "));
			    String curSeq = null;
			    //System.err.println(curName);
	    
			    if(len < 14){
				//second block contains the start of exon2
				if(len == 1){
				    curSeq = curline.substring(curline.indexOf("  ")+27).replaceAll("\\s+", "");
				    //System.err.println(curline + "\n" + curSeq);
				}else if(len == 13){ // 14th block contains the end of exon3
				    curSeq = curline.substring(curline.indexOf("  ")).substring(0, 15).replaceAll("\\s+", "");
				    //System.err.println(curline + "\n" + curSeq);
				    //System.err.println(curline + "\n" + curSeq + "|" + curSeq.length());
				}
				else
				    curSeq = curline.substring(curline.indexOf("  ")).replaceAll("\\s+","");
				if(hlas.get(curName) == null){
				    if(len != 1)
					System.err.println("Something is not right");
				    hlas.put(curName, new Sequence(curName, curSeq));
				    keys.add(curName);
				}else
				    hlas.get(curName).append(curSeq);
			    }
			}
		    }
		}else{//header section
		    //System.err.println(curline);
		    if(curline.equals("//")){
			br.readLine();
			inSequence = true;
			//System.err.println("HERE" + curline);
		    }
		}
	    }
	    br.close();
	}catch(IOException ioe){
	    ioe.printStackTrace();
	}
	//this.output();
    }

    public void processDQA1(String f){
	BufferedReader br = null;
	try{
	    br = new BufferedReader(new FileReader(f));
	    String curline = null;
	    boolean inSequence = false;
	    int len = 0;
	    while((curline=br.readLine())!=null){
		curline = curline.trim();
		/*if(curline.contains("A*80:03")){
		    System.err.println(curline);
		    String tmp = br.readLine();//System.err.println(br.readLine());
		    if(tmp.equals("//"))
			System.err.println("equals");
		    else
			System.err.println("Not euqls:" + tmp);
			}*/
		if(inSequence){//if we are in MSA
		    //System.err.println("inSequence");
		    if(curline.equals(""))// blank line means an end of a block
			len++;
		    else{//in MSA block
			//System.err.println("len\t-->\t" + len);
			if(len > 0){//first block we don't include since we are only targetting for exon2 and exon3
			    String curName = curline.substring(0,curline.indexOf("  "));
			    String curSeq = null;
			    //System.err.println(curName);
	    
			    if(len < 7){
				//second block contains the start of exon2
				if(len == 1){
				    curSeq = curline.substring(curline.indexOf("  ")+37).replaceAll("\\s+", "");
				    System.err.println(curline + "\n" + curSeq);
				}else if(len == 6){ // 15th block contains the end of exon3
				    curSeq = curline.substring(curline.indexOf("  ")).substring(0, 36).replaceAll("\\s+", "");
				    System.err.println(curline + "\n" + curSeq + "|" + curSeq.length());
				}
				else
				    curSeq = curline.substring(curline.indexOf("  ")).replaceAll("\\s+","");
				if(hlas.get(curName) == null){
				    if(len != 1)
					System.err.println("Something is not right");
				    hlas.put(curName, new Sequence(curName, curSeq));
				    keys.add(curName);
				}else
				    hlas.get(curName).append(curSeq);
			    }
			}
		    }
		}else{//header section
		    //System.err.println(curline);
		    if(curline.equals("//")){
			br.readLine();
			inSequence = true;
			//System.err.println("HERE" + curline);
		    }
		}
	    }
	    br.close();
	}catch(IOException ioe){
	    ioe.printStackTrace();
	}
	this.output();
    }



    private ArrayList<String> keys;
    private HashMap<String, Sequence> hlas;
}

class Sequence{

    private String id;
    private StringBuffer sequence;
    
    public Sequence(String id, String str){
	this.id = id;
	this.sequence = new StringBuffer(str.replaceAll("\\.+", ""));
    }
    
    public void append(String str){
	this.sequence.append(str.replaceAll("\\.+", ""));
    }
    
    public String getFastaStr(){
	return ">" + id + "\n"+ this.formatSequence();
    }
    private String formatSequence(){
	StringBuffer buffer = new StringBuffer();
	int lowInd = 0;
        for(int i=0; i<sequence.length();i++){
            if((i % 70) == 69){
                buffer.append(sequence.substring(lowInd,i+1));
                buffer.append("\n");
                lowInd += 70;
            }
        }
        if( (sequence.length()%70) != 0){
            buffer.append(sequence.substring(lowInd));
            buffer.append("\n");
        }
        return buffer.toString();
    }

}
