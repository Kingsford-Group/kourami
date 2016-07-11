import java.util.*;

public class NomG{
    
    private ArrayList<Group> groups;
    private HashMap<String, Group> allele2Group;
    
    public NomG(){
	this.groups = new ArrayList<Group>();
	this.allele2Group = new HashMap<String, Group>();
    }
    
}

class Group{

    private String hlaGeneName;

    private String groupname;
    
    private HashSet<String> set;
    
    public Group(String line){
	this.set = new HashSet<String>();
	this.process(line);
    }
    
    public void process(String line){
	String[] tokens = line.split(";");
	this.hlaGeneName = tokens[0].substring(0,tokens[0].indexOf("*"));
	if(tokens[2].trim().equals(""))
	    this.groupname = tokens[1];
	String[] elements = tokens[1].split("\\");
	for(String e : elements)
	    this.set.add(e);
    }
    
}

