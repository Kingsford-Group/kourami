import java.util.*;
import java.io.*;

public class NomG{
    
    private HashMap<String, ArrayList<Group>> hlagene2groups;

    //private ArrayList<Group> groups;
    private HashMap<String, Group> allele2Group;
    
    public NomG(){
	this.hlagene2groups = new HashMap<String, ArrayList<Group>>();
	//this.groups = new ArrayList<Group>();
	this.allele2Group = new HashMap<String, Group>();
    }

    
    //returns the list of groups belonging to HLA gene
    public ArrayList<Group> getGroups(String hgn){
	return this.hlagene2groups.get(hgn);
    }


    public void loadHlaGene2Groups(String nomGFile){
	BufferedReader br = null;
	String curline = null;
	try{
	    br = new BufferedReader(new FileReader(nomGFile));
	    while( (curline=br.readLine()) != null){
		if(curline.charAt(0) != '#'){
		    Group curgrp = new Group(curline, this);
		    String curhlagn = curgrp.getHLAGeneName();
		    ArrayList<Group> groups = this.hlagene2groups.get(curhlagn);
		    if(groups == null){
			groups = new ArrayList<Group>();
			this.hlagene2groups.put(curhlagn, groups);
		    }
		    groups.add(new Group(curline, this));
		}
	    }
	    br.close();
	}catch(IOException ioe){
	    ioe.printStackTrace();
	}
    }
    
    public void addToAllele2Group(String allele, Group g){
	this.allele2Group.put(allele, g);
    }
    
}

/*
 * G Group
 */
class Group{

    private String hlaGeneName;

    private String groupname;
    
    private HashSet<String> set;

    public String getGroupString(){
	return hlaGeneName + "*" + groupname;
    }

    public String getHLAGeneName(){
	return hlaGeneName;
    }
    
    public String getGroupName(){
	return this.groupname;
    }

    public String getFirstAllele(){
	return set.iterator().next();
    }
    
    public Group(String line, NomG nomG){
	this.set = new HashSet<String>();
	this.process(line, nomG);
    }
    
    public void process(String line, NomG nomG){
	String[] tokens = line.split(";");
	this.hlaGeneName = tokens[0].substring(0,tokens[0].indexOf("*"));
	if(tokens[2].trim().equals(""))
	    this.groupname = tokens[1];
	String[] elements = tokens[1].split("\\");
	for(String e : elements){
	    this.set.add(e);
	    nomG.addToAllele2Group(e, this);
	}
    }
    
}

