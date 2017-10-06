/*
Part of Kourami HLA typer/assembler
(c) 2017 by  Heewook Lee, Carl Kingsford, and Carnegie Mellon University.
See LICENSE for licensing.
*/
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceRecord;

import java.io.*;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.ArrayList;
import java.util.List;
import java.util.zip.GZIPInputStream;

import it.unimi.dsi.fastutil.objects.Object2IntOpenHashMap;

import org.apache.commons.cli.*;

public class HLA{

    public static boolean PRINT_G_GROUP_DB = false;
    /* graph mod stats */
    public static int NEW_NODE_ADDED = 0;
    public static int HOPPING = 0;
    public static int INSERTION_NODE_ADDED = 0;
    public static int INSERTION_WITH_NO_NEW_NODE = 0;
    public static int INSERTION = 0;
    /* end of graph mod stats */
    public static int READ_LENGTH = 100; //automatically gets set.
    public static double X_FACTOR = 4.0d/3.0d; //xFactor == 1 (a=4b), 4/3 (a=3b), 2 (a=2b)
    public static int SCORING_SCHEME = 4;//4 for APCUM+ISB

    public static LogHandler log;
    public static boolean DEBUG = false;
    public static boolean DEBUG3 = false;
    //output merged database MSA
    public static boolean OUTPUT_MERGED_MSA = false;

    //-o option
    public static String OUTPREFIX; // used for outfile names

    //-a option
    public static boolean TYPEADDITIONAL;

    //-d option
    public static String MSAFILELOC;
    public static String VERSION = "0.9.4";
    

    public HLA(String[] hlaList, String nomGFile){
	this.hlaName2Graph = new HashMap<String, HLAGraph>();
	this.hlaName2typingSequences = new HashMap<String, ArrayList<HLASequence>>();
	this.loadGraphs(hlaList, nomGFile);
    }

    //loads HLAGraphs as well as nomG typing sequences
    private void loadGraphs(String[] hlaList, String nomGFile){

	HLA.log.appendln("Merging HLA sequences and building HLA graphs");

	int i;
	NomG nomG = new NomG();
	nomG.loadHlaGene2Groups(nomGFile);
	
	String tmpDir = null;

	tmpDir = HLA.MSAFILELOC;
	
	for(i=0; i<hlaList.length; i++){
	    HLA.log.appendln("processing HLA gene:\t" + hlaList[i]);
	    //System.err.println("processing HLA gene:\t" + hlaList[i]);
	    MergeMSFs mm = new MergeMSFs();
	    if(!mm.merge(tmpDir + File.separator +  hlaList[i] + "_nuc.txt", tmpDir + File.separator + hlaList[i] + "_gen.txt", HLA.OUTPUT_MERGED_MSA)){
		HLA.log.appendln("ERROR in MSA merging. CANNOT proceed further. Exiting..");
		HLA.log.outToFile();
		System.exit(-1);
	    }
	    
	    this.hlaName2Graph.put(hlaList[i], new HLAGraph(mm.getListOfSequences(), hlaList[i]));
	    ArrayList<Group> groups = nomG.getGroups(hlaList[i]);
	    if(groups != null)
		this.hlaName2typingSequences.put(hlaList[i], mm.formDataBase(nomG.getGroups(hlaList[i])));
	    else
		this.hlaName2typingSequences.put(hlaList[i], mm.formDataBaseAll());
	    this.hlaName2Graph.get(hlaList[i]).setTypingSequences(this.hlaName2typingSequences.get(hlaList[i]));
	    if(HLA.OUTPUT_MERGED_MSA)
		this.outputTypingSequences(hlaList[i]);
	}
	HLA.log.appendln("Done building\t" + i + "\tgraphs.");
    }
    


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
    
    //kourami bam checker added
    private boolean checkHeader(SAMFileHeader header){
	List<SAMSequenceRecord> sequences = header.getSequenceDictionary().getSequences();
	HashSet<String> map = new HashSet<String>();
	
	//load kourami panel sequence names
	BufferedReader br;
	try{
	    br = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(HLA.MSAFILELOC + File.separator + "All_FINAL_with_Decoy.fa.gz"))));
	    String curline = "";
	    while((curline = br.readLine())!=null){
		if(curline.charAt(0) == ('>'))
		    map.add(curline.substring(1));
	    }
	    br.close();
	}catch(IOException ioe){
	    ioe.printStackTrace();
	}
	
	//check if input bam has sequences to kourami panel
	for(SAMSequenceRecord ssr : sequences){
	    if(!map.contains(ssr.getSequenceName()))
		return false;
	}
	return true;
    }


    public void loadReads(File[] bams) throws IOException{
	
	int count = 0;
	int numOp = 0;
	
	for(File bam : bams){
	    HLA.log.appendln("Loading reads from:\t" + bam.getName());
	    Object2IntOpenHashMap<String> readLoadingSet = new Object2IntOpenHashMap<String>();
	    readLoadingSet.defaultReturnValue(0);
	    
	    final SamReader reader = SamReaderFactory.makeDefault().open(bam);
	    
	    //Kourami bam checker added
	    if(!checkHeader(reader.getFileHeader())){
		HLA.log.appendln("Unexpected BAM :\t"+ bam.getName() 
				+"\nThe input BAM MUST be aligned to the set of IMGT/HLA alleles in " + HLA.MSAFILELOC + "\n" 
				+ "Please use the recommended preprocessing steps explained on the github page:\n"
				+ "https://github.com/Kingsford-Group/kourami");
		System.err.println("Unexpected BAM :\t"+ bam.getName() 
				   +"\nThe input BAM MUST be aligned to the set of IMGT/HLA alleles in " + HLA.MSAFILELOC + "\n" 
				   + "Please use the recommended preprocessing steps explained on the github page:\n"
				   + "https://github.com/Kingsford-Group/kourami");
		HLA.log.outToFile();
		System.exit(1);
	    }

	    for(final SAMRecord samRecord : reader){
		if(count == 0){
		    HLA.READ_LENGTH = samRecord.getReadLength();
		    HLA.log.appendln("Setting HLA.READ_LEGNTH = " + HLA.READ_LENGTH);
		}
		//added checking to process reads matching to HLA-type sequences
		//discarding decoy hits (DQB2, DQA2)
		boolean qc = false;
		if( (samRecord.getReferenceName().indexOf("*") > -1) 
		    && !samRecord.getReadUnmappedFlag() 
		    && !samRecord.isSecondaryOrSupplementary() 
		    && !this.startWIns(samRecord)){
		    count++;
		    if(samRecord.getReadPairedFlag())
			numOp += processRecord(samRecord, readLoadingSet);
		    else
			numOp += processRecordUnpaired(samRecord);
		}
		
		if(HLA.DEBUG && count%10000 == 0)
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
		
		readLoadingSet.put(sr.getReadName(), HLA.readNum);
		HLA.readNum++;
	    }else
		readnum = sr.getFirstOfPairFlag() ? readnum : 0-readnum;

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
	boolean readdebug = false;
	if(readdebug){
	    HLA.log.appendln(sr.getSAMString());
	    HLA.log.appendln("EffectiveLen:\t" + effectiveLen);
	    HLA.log.appendln("ReadLen:\t" + rLen);
	}
	Integer i = sr.getIntegerAttribute("NM");
	int nm = 0;
	if(i!=null)
	    nm = i.intValue();
	if(readdebug)
	    HLA.log.appendln("NM=\t" + nm);
	if(nm < 16){
	    if(readdebug)
		HLA.log.appendln("PASSWED QC");
	    return true;
	}
	if(readdebug){
	    HLA.log.appendln("FAILED QC");
	    HLA.log.appendln(sr.getSAMString());
	}
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

    public void printWeights(String[] list){
	for(String g:list)
	    this.hlaName2Graph.get(g).traverseAndWeights();
    }
    
    public void printBoundaries(String[] list){
	for(String g:list)
	    this.hlaName2Graph.get(g).getRefAllele().printBoundaries();
    }

    public void removeUnused(String[] list){
	for(String g:list)
	    this.hlaName2Graph.get(g).removeUnused();
    }

    public void flattenInsertionNodes(String[] list){
	for(String g:list)
	    this.hlaName2Graph.get(g).flattenInsertionNodes();
    }

    public void printStartEndNodes(String[] list){
	for(String g:list)
	    this.hlaName2Graph.get(g).printStartEndNodeInfo();
    }
    
    public void countBubbles(String[] list){
	for(String g:list)
	    this.hlaName2Graph.get(g).countBubbles();
    }
    
    public void countBubblesAndMerge(String[] list, StringBuffer rb){
	for(String g:list)
	    this.hlaName2Graph.get(g).countBubblesAndMerge(rb);
    }

    public void countStems(String[] list){
	for(String g:list)
	    this.hlaName2Graph.get(g).countStems();
    }
    
    public void removeStems(String[] list){
	for(String g:list)
	    this.hlaName2Graph.get(g).removeStems();
    }

    public void writeResults(StringBuffer rb, BufferedWriter resultWriter){
	try{
	    resultWriter.write(rb.toString());
	    resultWriter.close();
	}catch(IOException ioe){
	    ioe.printStackTrace();
	}
    }
    
    private static Options createHelpOption(){
	Options options = new Options();
	Option help = new Option("h", "help", false, "print this message");
	options.addOption(help);
	return options;
    }

    private static Options createOption(){
	Options options = new Options();
	Option help = new Option("h", "help", false, "print this message");

	Option buildFromMSA = Option.builder("d")
	    .longOpt("msaDirectory")
	    .required(true)
	    .argName("path")
	    .hasArg()
	    .desc("build HLA-Graph from gen and nuc MSAs provided by IMGT/HLA DB from given directory (required)")
	    .build();
	
	Option outfile = Option.builder("o")
	    .longOpt("outfilePrefix")
	    .required(true)
	    .hasArg()
	    .desc("use given outfile prefix for all output files (required)")
	    .argName("outfile")
	    .build();
	
	Option additionalLoci = Option.builder("a")
	    .longOpt("additionalLoci")
	    .required(false)
	    .hasArg(false)
	    .desc("type additional loci (optional)")
	    .build();

	//options.addOption(help);
	options.addOption(buildFromMSA);
	options.addOption(outfile);
	options.addOption(additionalLoci);
	
	return options;
    }

    private static void help(Options options){
	//String R = "\u001B[30m";
	HelpFormatter formatter = new HelpFormatter();
	formatter.setDescPadding(0);
	String header = "\n"
	    + "Program: Kourami - Graph-guided assembly of HLA typing exons\n"
	    + "Version: " + HLA.VERSION + "\n"
	    + "Contact: Heewook Lee <heewookl@cs.cmu.edu>\n\n"
	    + "Usage: java -jar <PATH_TO>/Kourami.jar [options] <bam-1> ... <bam-n>\n\n"
	    + "   -h,--help                      print this message\n";
	
	String footer = "\n";
	System.err.println(header);
	PrintWriter tmp = new PrintWriter(System.err);
	formatter.printOptions(tmp, 80, options, 3, 3);
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
	Options helponlyOpts = HLA.createHelpOption();
	String[] bams = null;
	CommandLine line = null;
	boolean exitRun = false;
	try{
	    CommandLine helpcheck = new DefaultParser().parse(helponlyOpts, args, true);
	    if(helpcheck.getOptions().length > 0)
		HLA.help(options);
	    else{
		line = parser.parse( options, args);
		if(line.hasOption("h"))//help"))
		    HLA.help(options);
		else{
		    if(line.hasOption("a"))
			HLA.TYPEADDITIONAL = true;

		    HLA.OUTPREFIX = line.getOptionValue("o");//outfilePrefix");
		    String tmploc = line.getOptionValue("d");//msaDirectory");
		    HLA.MSAFILELOC = tmploc;
		    if(tmploc.endsWith(File.separator))
			HLA.MSAFILELOC = tmploc.substring(0,tmploc.length()-1);
		    if(! new File(HLA.MSAFILELOC).exists() || ! new File(HLA.MSAFILELOC).isDirectory()){
			System.err.println("Given msaDirectory: " + HLA.MSAFILELOC + "\t does NOT exist or is NOT a directory.");
			exitRun = true;
		    }else if(! new File(HLA.MSAFILELOC + File.separator + "hla_nom_g.txt").exists()){
			System.err.println("hla_nom_g.txt NOT FOUND in " + HLA.MSAFILELOC );
			System.err.println("Please download hla_nom_g.txt from the same IMGT Release as msa files.");
			exitRun = true;
		    }
		}
		bams = line.getArgs();
		if(bams[bams.length - 1].equals("DEBUG1228")){
		    String[] tmpbams = new String[bams.length - 1]; 
		    for(int i=0;i<bams.length-1;i++)
			tmpbams[i] = bams[i];
		    bams = tmpbams;
		    HLA.DEBUG = true;
		}else
		    HLA.DEBUG = false;
		if(bams.length <1)
		    throw new ParseException("At least 1 bam file is required. See Usage:");
		else{
		    for(String b : bams)
			if(! new File(b).exists()){
			    System.err.println("Input bam : " + b + " DOES NOT exist. Please check the bam exists.");
			    exitRun = true;
			}
		}
		    
	    }
	    if(exitRun)
		throw new ParseException("Exitting . . .");
	}catch(ParseException e){
	    System.err.println(e.getMessage());
	    //System.err.println("Failed to parse command line args. Check usage.");
	    HLA.help(options);
	}
	
	String[] list = {"A" , "B" , "C" , "DQA1" , "DQB1" , "DRB1"};

	String[] extList = {"A" , "B" , "C" , "DQA1" , "DQB1" , "DRB1", "DOA", "DMA", "DMB"
			    ,"DPA1", "DPB1", "DRA",  "F", "G" , "H", "J" , "L"};
	//,"DPA1", "DPB1", "DRA",  "DRB4", "F", "G" , "H", "J" ,"K", "L", "V"};
	//,"DPA1", "DPB1", "DRA", "DRB3", "DRB4", "F", "G" , "H", "J" ,"K", "L", "V"};
	
	if(HLA.TYPEADDITIONAL)
	    list = extList;

	File[] bamfiles = new File[bams.length];

	for(int i=0;i<bams.length; i++)
	    bamfiles[i] = new File(bams[i]);
	
	//check if <HLA.OUTPREFIX>.result is writable
	//if not exit.
	BufferedWriter resultWriter = null;
	try{
	    resultWriter = new BufferedWriter(new FileWriter(HLA.OUTPREFIX + ".result"));
	}catch(IOException ioe){
	    ioe.printStackTrace();
	    System.err.println("\n\n>>> CANNOT open output file: " + HLA.OUTPREFIX + ".result <<<\n\n");
	    HLA.help(options);
	}

	
	HLA.log = new LogHandler();
	for(int i =0; i<args.length;i++)
	    HLA.log.append(" " + args[i]);
	HLA.log.appendln();

	try{
	    System.err.println("----------------REF GRAPH CONSTRUCTION--------------");
	    HLA.log.appendln("----------------REF GRAPH CONSTRUCTION--------------");
	    HLA hla = new HLA(list, HLA.MSAFILELOC + File.separator + "hla_nom_g.txt");

	    //1. bubble counting before loading reads.
	    //System.err.println("----------------BUBBLE COUNTING: REF GRAPH--------------");
	    //HLA.log.appendln("----------------BUBBLE COUNTING: REF GRAPH--------------");
	    //hla.countStems();
	    
	    System.err.println("---------------- READ LOADING --------------");
	    HLA.log.appendln("---------------- READ LOADING --------------");
	    
	    hla.loadReads(bamfiles); 
	    
	    System.err.println("---------------- GRAPH CLEANING --------------");
	    HLA.log.appendln("---------------- GRAPH CLEANING --------------");
	    	    	    
	    hla.flattenInsertionNodes(list);
	    hla.removeUnused(list);
	    hla.removeStems(list);
	    
	    /*updating error prob*/
	    hla.updateErrorProb();
	    
	    hla.log.flush();
	    
	    StringBuffer resultBuffer = new StringBuffer();
	    
	    HLA.DEBUG3 = HLA.DEBUG;

	    hla.countBubblesAndMerge(list, resultBuffer);
	    
	    hla.writeResults(resultBuffer, resultWriter);
	}catch(Exception e){
	    e.printStackTrace();
	    HLA.log.outToFile();
	    System.exit(-1);
	}
	/*printingWeights*/
	//hla.printWeights();
	HLA.log.outToFile();
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
}
