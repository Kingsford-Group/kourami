public class GetRanking{
    private HashMap<String,String> allele2G_A;
    private HashMap<String,String> allele2G_B;
    private HashMap<String,String> allele2G_C;
    private HashMap<String,String> allele2G_DQA1;
    private HashMap<String,String> allele2G_DQB1;
    private HashMap<String,String> allele2G_DRB1;


    public GetRanking(){
	this.allele2G_A = new HashMap<String, String>();
	this.allele2G_B = new HashMap<String, String>();
	this.allele2G_C = new HashMap<String, String>();
	this.allele2G_DQA1 = new HashMap<String, String>();
	this.allele2G_DQB1 = new HashMap<String, String>();
	this.allele2G_DRB1 = new HashMap<String, String>();
    }

    public void loadHash(String nomgfile){
	BufferedReader br = null;
	try{
	    br = new BufferedReader(new FileReader(nomgfile));
	    String curline = "";
	    while((curline=br.readLine())!=null){
		if(curline.startsWith("#"))
		    ;
		else{
		    String[] tokens = curline.split(";");
		    if(tokens[0].equals("A*")){
			this.add2Hash(allele2G_A, tokens[1], tokens[2], tokens[0]);
		    }else if(tokens[0].equals("B*")){
			this.add2Hash(allele2G_B, tokens[1], tokens[2], tokens[0]);
		    }else if(tokens[0].equals("C*")){
			this.add2Hash(allele2G_C, tokens[1], tokens[2], tokens[0]);
		    }else if(tokens[0].equals("DQA1*")){
			this.add2Hash(allele2G_DQA1, tokens[1], tokens[2], tokens[0]);
		    }else if(tokens[0].equals("DQB1*")){
			this.add2Hash(allele2G_DQB1, tokens[1], tokens[2], tokens[0]);
		    }else if(tokens[0].equals("DRB1*")){
			this.add2Hash(allele2G_DRB1, tokens[1], tokens[2], tokens[0]);
		    }
		}
	    }
	    br.close();
	}catch(IOException ioe){
	    ioe.printStackTrace();
	}
    }
    
    private void add2Hash(HashMap<String,String> hash, String alleles, String g, String gene){
	if(g.length() > 0){//redundant alleles
	    String[] tokens = alleles.split("/");
	    for(int i=0;i<tokens.length;i++)
		hash.put(gene+tokens[i], gene+g);
	    
	}else{
	    hash.put(gene+alleles, gene+alleles);
	}
    }
    
}
