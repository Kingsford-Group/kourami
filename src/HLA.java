import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

import java.io.*;
import java.util.HashMap;
import java.util.Iterator;
import java.util.ArrayList;

import it.unimi.dsi.fastutil.objects.Object2IntOpenHashMap;

import org.apache.commons.cli.*;

public class HLA{

    public static boolean PRINT_G_GROUP_DB = false;
    public static boolean ALIGNMATCHING_WHEN_TIED = false;
    public static int NEW_NODE_ADDED = 0;
    public static int HOPPING = 0;
    public static int INSERTION_NODE_ADDED = 0;
    public static int INSERTION_WITH_NO_NEW_NODE = 0;
    public static int INSERTION = 0;
    public static int READ_LENGTH = 100; //automatically gets set.
    public static double X_FACTOR = 4.0d/3.0d; //xFactor == 1 (a=4b), 4/3 (a=3b), 2 (a=2b)
    
    public static LogHandler log;
    public static boolean DEBUG = true;

    public static boolean OUTPUT_MERGED_MSA = false;

    public static String OUTPREFIX; // used for outfile names
    //public static boolean SERIALIZEIT = false;
    //public static String SERIALIZEFILE;
    public static String MSAFILELOC;
    //public static String PREBUILTFILE = ".." + File.separator + "db" + File.separator + "3240.db";
    public static String VERSION = "0.1.0";
    

    public HLA(String[] hlaList, String nomGFile){
	this.hlaName2Graph = new HashMap<String, HLAGraph>();
	this.hlaName2typingSequences = new HashMap<String, ArrayList<HLASequence>>();
	this.loadGraphs(hlaList, nomGFile);
    }

        //loads HLAGraphs as well as nomG typing sequences
    private void loadGraphs(String[] hlaList, String nomGFile){
	//if(HLA.PREBUILTFILE == null)
	HLA.log.appendln("Merging HLA sequences and building HLA graphs");
	//else
	    //	    HLA.log.appendln("Loading prebuilt MSAs and building HLA graphs");
	int i;
	NomG nomG = new NomG();
	nomG.loadHlaGene2Groups(nomGFile);
	
	String tmpDir = null;
	//if(HLA.PREBUILTFILE == null)
	tmpDir = HLA.MSAFILELOC;//"../db";
	
	for(i=0; i<hlaList.length; i++){
	    HLA.log.appendln("processing HLA gene:\t" + hlaList[i]);
	    MergeMSFs mm = new MergeMSFs();
	    //if(HLA.PREBUILTFILE == null){
	    //mm = new MergeMSFs();
	    if(!mm.merge(tmpDir + File.separator +  hlaList[i] + "_nuc.txt", tmpDir + File.separator + hlaList[i] + "_gen.txt", HLA.OUTPUT_MERGED_MSA)){
		System.err.println("ERROR in MSA merging. CANNOT proceed further. Exiting..");
		System.exit(-1);
	    }
		//if(HLA.SERIALIZEIT)
		//  this.serializeMergeMSFs(mm, HLA.SERIALIZEFILE + "." + hlaList[i]);
		
	    //}else
	    //	mm = this.deserializeMergeMSFs(HLA.PREBUILTFILE + "." + hlaList[i]);
	    
	    this.hlaName2Graph.put(hlaList[i], new HLAGraph(mm.getListOfSequences(), hlaList[i]));
	    this.hlaName2typingSequences.put(hlaList[i], mm.formDataBase(nomG.getGroups(hlaList[i])));
	    this.hlaName2Graph.get(hlaList[i]).setTypingSequences(this.hlaName2typingSequences.get(hlaList[i]));
	    if(HLA.OUTPUT_MERGED_MSA)
		this.outputTypingSequences(hlaList[i]);
	}
	HLA.log.appendln("Done building\t" + i + "\tgraphs.");
    }
    

    /*
    private void loadGraphs(String[] hlaList, String nomGFile){
	String tmpDir = "../db";
	HLA.log.appendln("Merging HLA sequences and building HLA graphs");
	int i;
	NomG nomG = new NomG();
	nomG.loadHlaGene2Groups(nomGFile);
	for(i=0; i<hlaList.length; i++){
	    HLA.log.appendln("processing HLA gene:\t" + hlaList[i]);
	    MergeMSFs mm = new MergeMSFs();
	    if(!mm.merge(tmpDir + File.separator +  hlaList[i] + "_nuc.txt", tmpDir + File.separator + hlaList[i] + "_gen.txt", HLA.OUTPUT_MERGED_MSA)){
		System.err.println("ERROR in MSA merging. CANNOT proceed further. Exiting..");
		System.exit(-1);
	    }
	    this.hlaName2Graph.put(hlaList[i], new HLAGraph(mm.getListOfSequences()));
	    this.hlaName2typingSequences.put(hlaList[i], mm.formDataBase(nomG.getGroups(hlaList[i])));
	    this.hlaName2Graph.get(hlaList[i]).setTypingSequences(this.hlaName2typingSequences.get(hlaList[i]));
	    this.outputTypingSequences(hlaList[i]);
	}
	HLA.log.appendln("Done building\t" + i + "\tgraphs.");
	}*/

    /*
    public void serializeMergeMSFs(MergeMSFs mm, String outFileName){
	try{
	    HLA.log.appendln(">>>>>>>>> Serializing MergeMSFs <<<<<<<<<<");
	    FileOutputStream f = new FileOutputStream(new File(outFileName));
	    ObjectOutputStream o = new ObjectOutputStream(f);
	    o.writeObject(mm);
	    o.close();
	    f.close();
	    
	}catch(FileNotFoundException e){
	    HLA.log.appendln("Error in MergeMSFs Serialization: File not found");
	    HLA.log.flush();
	    e.printStackTrace();
	}catch(IOException e){
	    HLA.log.appendln("Error in MergeMSFs Serialization: error in stream initialization");
	    HLA.log.flush();
	    e.printStackTrace();
	}//catch(ClassNotFoundException e){
	 //   e.printStackTrace();
	//}
    }
    
    public MergeMSFs deserializeMergeMSFs(String serializedF){
	MergeMSFs mm = null;
	try{
	    HLA.log.appendln("<<<<<<<<< De-Serializaing MErgeMSFs >>>>>>>>>>>");
	    FileInputStream fi = new FileInputStream(new File(serializedF));
	    ObjectInputStream oi = new ObjectInputStream(fi);
	    mm = (MergeMSFs) oi.readObject();
	    oi.close();
	    fi.close();
	    return mm;
	}catch(FileNotFoundException e){
	    HLA.log.appendln("Error in MergeMSFs De-serialization: File not found");
	    HLA.log.flush();
	}catch(IOException e){
	    HLA.log.appendln("Error in MergeMSFs De-serialization: error in stream initialization");
	    HLA.log.flush();
	}catch(ClassNotFoundException e){
	  e.printStackTrace();
	}
	return null;
    }
    */
    public void outputTypingSequences(String hgn){
	ArrayList<HLASequence> typingSeqs = this.hlaName2typingSequences.get(hgn);

	BufferedWriter bw = null;
	try{
	    bw = new BufferedWriter(new FileWriter(hgn + "_typingDB.fa"));
	    for(HLASequence h : typingSeqs){
		bw.write(Bubble.stripPadding(h.toString()));
		//bw.write(">" + h.getGroup().getGroupString() + "\n");
		//bw.write();
	    }
	    bw.close();
	}catch(IOException ioe){
	    ioe.printStackTrace();
	}
    }

    public void loadReads(File[] bams) throws IOException{
	
	int count = 0;
	int numOp = 0;
	
	for(File bam : bams){
	    HLA.log.appendln("Loading reads from:\t" + bam.getName());
	    Object2IntOpenHashMap<String> readLoadingSet = new Object2IntOpenHashMap<String>();
	    //if(pairedend){
	    //readLoadingSet = new Object2IntOpenHashMap<String>();
	    readLoadingSet.defaultReturnValue(0);
	    //}
	    
	    final SamReader reader = SamReaderFactory.makeDefault().open(bam);
	    for(final SAMRecord samRecord : reader){
		if(count == 0){
		    HLA.READ_LENGTH = samRecord.getReadLength();
		    HLA.log.appendln("Setting HLA.READ_LEGNTH = " + HLA.READ_LENGTH);
		}
		    //HLA.READ_LENGTH = (samRecord.getReadLength() > 300) ? 100 : samRecord.getReadLength();
		//System.out.println(samRecord.getCigarString());
		//samRecord
		
		//added checking to process reads matching to HLA-type sequences
		//discarding decoy hits (DQB2, DQA2)
		boolean qc = false;
		if( (samRecord.getReferenceName().indexOf("*") > -1) 
		    && !samRecord.getReadUnmappedFlag() 
		    && !samRecord.isSecondaryOrSupplementary() 
		    && !this.startWIns(samRecord)){
		    //		    && (qc==this.qcCheck(samRecord)) ){
		    count++;
		    if(samRecord.getReadPairedFlag())
			numOp += processRecord(samRecord, readLoadingSet);
		    else
			numOp += processRecordUnpaired(samRecord);
		}/*else{
		    if(!qc){
			HLA.log.appendln("SKIPPING LOW-QUAL READS");
			HLA.log.appendln(samRecord.getSAMString());
		    }
		    }*/
		if(count%10000 == 0)
		    HLA.log.appendln("Processed 10000 reads...");
	    }
	    reader.close();
	}
	HLA.log.appendln("Loaded a total of " + count + " mapped reads.");
	HLA.log.appendln("A total of " + numOp + " bases");
    }
    
    public void updateErrorProb(){
	HLA.log.appendln("------------ UPDATING error probabilities of each edge ---------");
	Iterator itr = this.hlaName2Graph.keySet().iterator();
	while(itr.hasNext()){
	    this.hlaName2Graph.get(itr.next()).updateEdgeWeightProb();
	}
	HLA.log.appendln("------------     DONE UPDATING error probabilities     ---------");
    }
    
    //assume interleaved SAMRecord
    public int processRecord(SAMRecord sr, Object2IntOpenHashMap<String> readLoadingSet){
	int totalOp = 0;
	String hlagene = HLA.extractHLAGeneName(sr.getReferenceName());
	HLAGraph hg = this.hlaName2Graph.get(hlagene);
	//hg.traverse();
	if(hg != null){
	    if(hg.isClassI()){
		boolean qc = this.qcCheck(sr);
		if(!qc)
		    return 0;
	    }
	    int readnum = readLoadingSet.getInt(sr.getReadName());
	    //no such read has been read. return value of 0 means the hashSet doesn't have the read
	    if(readnum == 0){
		readnum = sr.getFirstOfPairFlag() ? HLA.readNum : 0-HLA.readNum;
		//(132)  (-5695)  (-5137)  (5137)  (-2567)  (-1363)  (5143)
		
		readLoadingSet.put(sr.getReadName(), HLA.readNum);
		HLA.readNum++;
	    }else
		readnum = sr.getFirstOfPairFlag() ? readnum : 0-readnum;
	    
	    /*if(readnum == 2197 || readnum == -2197 || readnum == 2199 || readnum == -2199 || readnum == 2196 || readnum == -2196 || readnum ==2198 || readnum ==-2198 || readnum ==267 || readnum == -267){
		HLA.log.appendln("readnumPROBLEM( " + readnum + "):" + sr.getReadName());
		HLA.log.appendln(sr.getSAMString());
		}*/
	    /*
	    if(readnum == 1365 || readnum == -3991 || readnum == 5164 || readnum == -415 || readnum == 780 || readnum == -631 ){//if(readnum == 5137 || readnum == -5137 || readnum == 5143 || readnum == 132 || readnum == -5695 || readnum == -2567 || readnum == -1363){
	    HLA.log.appendln("readnumPROBLEM( " + readnum + "):" + sr.getReadName());
		}*/

	    totalOp += hg.addWeight(sr, readnum);//HLA.readNum);
	    //HLA.readNum++;
	}else{
	    ;//HLA.log.appendln("UNKNOWN HLA GENE: " + hlagene);
	}
	return totalOp;
    }
    
    public boolean startWIns(SAMRecord sr){
	Cigar cigar = sr.getCigar();
	if(cigar == null){
	    return true;
	}else{
	    CigarOperator op = cigar.getCigarElements().get(0).getOperator();
	    if(op == CigarOperator.I){
		if(HLA.DEBUG)
		    HLA.log.appendln("SKIPPING(Start with Insertion):\t" + sr.getReadName());
		return true;
		
	    }
	}
	return false;
    }

    public boolean qcCheck(SAMRecord sr){
	boolean debug = HLA.DEBUG;
	Cigar cigar = sr.getCigar();
	int rLen = sr.getReadLength();
	int effectiveLen = 0;
	if(cigar==null) 
	    return false;
	else{
	    for(final CigarElement ce : cigar.getCigarElements()){
		CigarOperator op = ce.getOperator();
		int cigarLen = ce.getLength();
		switch(op)
		    {
		    case M:
			{
			    effectiveLen += cigarLen;
			    break;
			}
		    case I:
			{
			    effectiveLen += cigarLen;
			    break;
			}
		    default: 
			break;
		    }
	    }
	}
	if(debug){
	    HLA.log.appendln(sr.getSAMString());
	    HLA.log.appendln("EffectiveLen:\t" + effectiveLen);
	    HLA.log.appendln("ReadLen:\t" + rLen);
	}
	Integer i = sr.getIntegerAttribute("NM");
	int nm = 0;
	if(i!=null)
	    nm = i.intValue();
	if(debug)
	    HLA.log.appendln("NM=\t" + nm);
	//if((effectiveLen*1.0d)/(rLen*1.0d) > 0.7d && (nm < 3 || ((nm*1.0d)/(effectiveLen*1.0d)) < 0.05d)){
	if(nm < 16){
	    if(debug)
		HLA.log.appendln("PASSWED QC");
	    return true;
	}
	if(debug)
	    HLA.log.appendln("FAILED QC");
	HLA.log.appendln(sr.getSAMString());
	return false;
    }


    public int processRecordUnpaired(SAMRecord sr){
	int totalOp = 0;
	String hlagene = HLA.extractHLAGeneName(sr.getReferenceName());
	HLAGraph hg = this.hlaName2Graph.get(hlagene);
	//hg.traverse();
	if(hg != null){
	    if(hg.isClassI()){
		boolean qc = this.qcCheck(sr);
		if(!qc)
		    return 0;
	    }
	    totalOp += hg.addWeight(sr, HLA.readNum);
	    HLA.readNum++;
	}
	return totalOp;
    }

    /*    
    public int processRecord(SAMRecord sr){
	int totalOp = 0;
	String hlagene = HLA.extractHLAGeneName(sr.getReferenceName());
	HLAGraph hg = this.hlaName2Graph.get(hlagene);
	//hg.traverse();
	if(hg != null){
	    totalOp += hg.addWeight(sr, HLA.readNum);
	    HLA.readNum++;
	}else{
	    ;
	}
	return totalOp;
    }
    */
    /*
    public void setNames(){
	this.hlaName2Graph.get("A").setHLAGeneName("A");
	this.hlaName2Graph.get("B").setHLAGeneName("B");
	this.hlaName2Graph.get("C").setHLAGeneName("C");
	this.hlaName2Graph.get("DQA1").setHLAGeneName("DQA1");
	this.hlaName2Graph.get("DQB1").setHLAGeneName("DQB1");
	this.hlaName2Graph.get("DRB1").setHLAGeneName("DRB1");
    }
    */
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

    public void removeUnused(){
	this.hlaName2Graph.get("A").removeUnused();
	this.hlaName2Graph.get("B").removeUnused();
	this.hlaName2Graph.get("C").removeUnused();
	this.hlaName2Graph.get("DQA1").removeUnused();
	this.hlaName2Graph.get("DQB1").removeUnused();
	this.hlaName2Graph.get("DRB1").removeUnused();
    }

    public void flattenInsertionNodes(){
	this.hlaName2Graph.get("A").flattenInsertionNodes();
	this.hlaName2Graph.get("B").flattenInsertionNodes();
	this.hlaName2Graph.get("C").flattenInsertionNodes();
	this.hlaName2Graph.get("DQA1").flattenInsertionNodes();
	this.hlaName2Graph.get("DQB1").flattenInsertionNodes();
	this.hlaName2Graph.get("DRB1").flattenInsertionNodes();
    }


    public void printStartEndNodes(){
	this.hlaName2Graph.get("A").printStartEndNodeInfo();
	this.hlaName2Graph.get("B").printStartEndNodeInfo();
	this.hlaName2Graph.get("C").printStartEndNodeInfo();
	this.hlaName2Graph.get("DQA1").printStartEndNodeInfo();
	this.hlaName2Graph.get("DQB1").printStartEndNodeInfo();
	this.hlaName2Graph.get("DRB1").printStartEndNodeInfo();
    }

    public void countBubbles(){
	this.hlaName2Graph.get("A").countBubbles();
	this.hlaName2Graph.get("B").countBubbles();
	this.hlaName2Graph.get("C").countBubbles();
	this.hlaName2Graph.get("DQA1").countBubbles();
	this.hlaName2Graph.get("DQB1").countBubbles();
	this.hlaName2Graph.get("DRB1").countBubbles();
    }

    public void countBubblesAndMerge(StringBuffer rb){
	this.hlaName2Graph.get("A").countBubblesAndMerge(rb);
	this.hlaName2Graph.get("B").countBubblesAndMerge(rb);
	this.hlaName2Graph.get("C").countBubblesAndMerge(rb);
	this.hlaName2Graph.get("DQA1").countBubblesAndMerge(rb);
	this.hlaName2Graph.get("DQB1").countBubblesAndMerge(rb);
	this.hlaName2Graph.get("DRB1").countBubblesAndMerge(rb);
    }
    /*
    public void countBubblesAndMerge(){
	this.hlaName2Graph.get("DQA1").countBubblesAndMerge();
	this.hlaName2Graph.get("DQB1").countBubblesAndMerge();
	this.hlaName2Graph.get("DRB1").countBubblesAndMerge();
	this.hlaName2Graph.get("A").countBubblesAndMerge();
	this.hlaName2Graph.get("B").countBubblesAndMerge();
	this.hlaName2Graph.get("C").countBubblesAndMerge();
	}*/
    
    public void countStems(){
	this.hlaName2Graph.get("A").countStems();
	this.hlaName2Graph.get("B").countStems();
	this.hlaName2Graph.get("C").countStems();
	this.hlaName2Graph.get("DQA1").countStems();
	this.hlaName2Graph.get("DQB1").countStems();
	this.hlaName2Graph.get("DRB1").countStems();
    }

    public void removeStems(){
	this.hlaName2Graph.get("A").removeStems();
	this.hlaName2Graph.get("B").removeStems();
	this.hlaName2Graph.get("C").removeStems();
	this.hlaName2Graph.get("DQA1").removeStems();
	this.hlaName2Graph.get("DQB1").removeStems();
	this.hlaName2Graph.get("DRB1").removeStems();
    }

    public void setFileName(String f){
	this.outfilename = f;
	this.hlaName2Graph.get("A").setFileName(f);
	this.hlaName2Graph.get("B").setFileName(f);
	this.hlaName2Graph.get("C").setFileName(f);
	this.hlaName2Graph.get("DQA1").setFileName(f);
	this.hlaName2Graph.get("DQB1").setFileName(f);
	this.hlaName2Graph.get("DRB1").setFileName(f);
    }


    public void writeResults(StringBuffer rb){
	BufferedWriter bw = null;
	try{
	    bw = new BufferedWriter(new FileWriter(this.outfilename + ".result"));
	    bw.write(rb.toString());
	    bw.close();
	}catch(IOException ioe){
	    ioe.printStackTrace();
	}
    }
    
    private static Options createOption(){
	Options options = new Options();
	
	Option help = new Option("help", "Print this message");
	//	Option buildFromMSA = new Option("buildFromMSA", "build HLAGraph from gen and nuc MSAs provided by IMGT/HLA DB");
	
	
	Option buildFromMSA = Option.builder("msaDirectory")
	    .required(true)
	    .argName("msaDirectory")
	    .hasArg()
	    .desc("build HLAGraph from gen and nuc MSAs provided by IMGT/HLA DB from given directory")
	    .build();
	//.create("db_filename");

	/*
	Option serialize = Option.builder( "serialize")
	    .hasArg()
	    .desc("serialize the constructed msa to given file")
	    .argName("file")
	    .build();
	
	Option usePrebuiltDB = Option.builder("usePrebuiltDB")
	    .hasArg()
	    .desc("use given prebuilt serialized DB (default: IMGT/HLA 3.24.0)")
	    .argName("usePrebuiltDB")
	    .build();
	*/
	
	Option outfile = Option.builder("outfilePrefix")
	    .required(true)
	    .hasArg()
	    .desc("use given outfile prefix for all output files")
	    .argName("outfile")
	    .build();
	
	options.addOption(help);
	options.addOption(buildFromMSA);
	//options.addOption(serialize);
	//options.addOption(usePrebuiltDB);
	options.addOption(outfile);
	
	return options;
    }

    private static void help(Options options){
	String R = "\u001B[30m";
	HelpFormatter formatter = new HelpFormatter();
	formatter.setDescPadding(0);
	String header = "\n"
	    + "Program: Kourami - Assembly of HLA typing exons\n"
	    + "Version: " + HLA.VERSION + "\n"
	    + "Contact: Heewook Lee <heewookl@cs.cmu.edu>\n\n"
	    + "Usage: java -jar <PATH_TO>/Kourami.jar [options] <bam-1> ... <bam-n>\n\n";
	
	String footer = "\n";
	//formatter.printHelp("Main", options);
	//formatter.printHelp(70, "K", header, options, footer, false);
	System.err.println(header);
	PrintWriter tmp = new PrintWriter(System.err);
	formatter.printOptions(tmp, 100, options, 3, 3);
	tmp.println("\n");
	tmp.println("            -hhy+.                o o       o o       o o o o       o o");
	tmp.println(".`           -syss:---.`        o     o o o     o o o         o o o     o o o");
	tmp.println(":+:`     .:/o+++++///ommy+`    o       _  __                               _");
	tmp.println("`yhs/..:osssooooo++++dmNNNdo`   o     | |/ /___  _   _ _ __ __ _ _ __ ___ (_)");
	tmp.println(" /syy///++++ooooooooodNMdNdmh: o      | ' // _ \\| | | | '__/ _` | '_ ` _ \\| |");
	tmp.println(" -do/` .://++++++++oodmmmmmmd-        | . \\ (_) | |_| | | | (_| | | | | | | |");
	tmp.println(" .+:     `.://///+///ommmmdy-         |_|\\_\\___/ \\__,_|_|  \\__,_|_| |_| |_|_|");
	tmp.println("  .          -syo----..``          ");
	tmp.println("            +y+.                \n\n");

	tmp.flush();
	tmp.close();
	System.exit(1);
    }

    public static void main(String[] args) throws IOException{
	
	CommandLineParser parser = new DefaultParser();
	Options options = HLA.createOption();
	String[] bams = null;
	CommandLine line = null;
	try{
	    line = parser.parse( options, args);
	    if(line.hasOption("help"))
		HLA.help(options);
	    else{
		HLA.OUTPREFIX = line.getOptionValue("outfilePrefix");
		String tmploc = line.getOptionValue("msaDirectory");
		HLA.MSAFILELOC = tmploc;
		if(tmploc.endsWith(File.separator))
		    HLA.MSAFILELOC = tmploc.substring(0,tmploc.length()-1);
		/*
		//can't have both options turned on
		if(line.hasOption("buildFromMSA") && line.hasOption("usePrebuiltDB")){
		    System.err.println("imcompatible options:");
		    HLA.help(options);
		}//else if(!line.hasOption("outfilePrefix"))
		//    System.err.println("-outfilePrefix <outfile> is required.");
		else{
		    HLA.OUTPREFIX = line.getOptionValue("outfilePrefix");
		    
		    //if building from gen and nuc files
		    if(line.hasOption("buildFromMSA")){
			HLA.PREBUILTFILE = null;
			String tmploc = line.getOptionValue("buildFromMSA");
			HLA.MSAFILELOC = tmploc;
			if(tmploc.endsWith(File.separator))
			    HLA.MSAFILELOC = tmploc.substring(0,tmploc.length()-1);
			if(line.hasOption("serialize")){
			    HLA.SERIALIZEIT = true;
			    HLA.SERIALIZEFILE = line.getOptionValue("serialize");
			}else
			    HLA.SERIALIZEIT = false;
		    }else if(line.hasOption("usePrebuiltDB"))
			HLA.PREBUILTFILE = line.getOptionValue("usePrebuiltDB");
			}*/
	    }
	    bams = line.getArgs();
	    if(bams.length <1)
		throw new ParseException("At least 1 bam file is required. See Usage:");
	}catch(ParseException e){
	    System.err.println(e.getMessage());
	    //System.err.println("Failed to parse command line args. Check usage.");
	    HLA.help(options);
	}
	
	String[] list = {"A" , "B" , "C" , "DQA1" , "DQB1" , "DRB1"};
	File[] bamfiles = new File[bams.length];

	for(int i=0;i<bams.length; i++)
	    bamfiles[i] = new File(bams[i]);
	String outfilename = HLA.OUTPREFIX;
	
	
	HLA.log = new LogHandler(outfilename);
	for(int i =0; i<args.length;i++)
	    HLA.log.append(" " + args[i]);
	HLA.log.appendln();

	try{
	    HLA hla = new HLA(list, HLA.MSAFILELOC + File.separator + "hla_nom_g.txt");
	    //sets HLA geneNames to each graph.
	    //hla.setNames();
	    //hla.printBoundaries();
	    //1. bubble counting before loading reads.
	    //System.err.println("----------------BUBBLE COUNTING: REF GRAPH--------------");
	    //HLA.log.appendln("----------------BUBBLE COUNTING: REF GRAPH--------------");
	    
	    //hla.countStems();
	    
	    System.err.println("---------------- READ LOADING --------------");
	    HLA.log.appendln("---------------- READ LOADING --------------");
	    
	    hla.loadReads(bamfiles); 
	    hla.setFileName(outfilename);
	    
	    System.err.println("---------------- GRAPH CLEANING --------------");
	    HLA.log.appendln("---------------- GRAPH CLEANING --------------");
	    	    	    
	    hla.flattenInsertionNodes();
	    hla.removeUnused();
	    hla.removeStems();
	    //hla.countStems();
	    
	    /*updating error prob*/
	    hla.updateErrorProb();
	    
	    hla.log.flush();
	    
	    StringBuffer resultBuffer = new StringBuffer();
	    
	    hla.countBubblesAndMerge(resultBuffer);
	    
	    hla.writeResults(resultBuffer);
	}catch(Exception e){
	    e.printStackTrace();
	    HLA.log.outToFile();
	    System.exit(-1);
	}
	/*printingWeights*/
	//hla.printWeights();
	HLA.log.outToFile();
	//public static int NEW_NODE_ADDED = 0;
	//public static int HOPPING = 0;
	//public static int INSERTION_NODE_ADDED = 0;
    	HLA.log.appendln("NEW_NODE_ADDED:\t" + HLA.NEW_NODE_ADDED);
	HLA.log.appendln("HOPPPING:\t" + HLA.HOPPING);
	HLA.log.appendln("INSERTION_NODE_ADDED:\t" + HLA.INSERTION_NODE_ADDED);
	HLA.log.appendln("INSERTION_WITH_NO_NEW_NODE:\t" + HLA.INSERTION_WITH_NO_NEW_NODE);
	HLA.log.appendln("INSERTION_COUNTS:\t" + HLA.INSERTION);
    }

    private static String extractHLAGeneName(String g){
	//if(g.indexOf("*") < 0)
	//return null;
	return g.substring(0,g.indexOf("*"));
    }
    
    
    public static int readNum = 1;
    private HashMap<String, HLAGraph> hlaName2Graph;
    private HashMap<String, ArrayList<HLASequence>> hlaName2typingSequences;
    private String outfilename;
}
