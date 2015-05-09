package jped;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;
import java.util.Collections;

/** Command line parser from apache */
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.GnuParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.commons.cli.Parser;

/** Classes from the Jannovar library. */
import jannovar.annotation.Annotation;
import jannovar.annotation.AnnotationList;
import jannovar.common.VariantType;
import jannovar.exception.JannovarException;
import jannovar.exome.Variant;
import jannovar.genotype.GenotypeCall;
import jannovar.io.PedFileParser;
import jannovar.io.SerializationManager;
import jannovar.io.VCFReader;
import jannovar.pedigree.Pedigree;
import jannovar.reference.Chromosome;
import jannovar.reference.TranscriptModel;

/**
 * Example Jannovar program designed to demonstrate how to integrate simple pedigree
 * filtering with exome variant annotation.
 * @author Peter Robinson
 * @version 1.2 (December 27, 2013)
 */
public class JPed {

    /** Map of Chromosomes (Jannovar object storing all transcripts from the Chromosome). */
    private HashMap<Byte,Chromosome> chromosomeMap=null;
     /** List of variants from input file to be analysed. */
    private ArrayList<Variant> variantList=null;
    /** List of all sample names of VCF file */
    private ArrayList<String> sampleNames = null;
    /** Pedigree of the persons whose samples were sequenced. Created on the basis of a PED file
	for multisample VCF files, or as a default single-sample Pedigree for single-sample VCF
	files. */
    private Pedigree pedigree = null;
    /** Map of all Genes encountered during parsing of VCF file. The key is the gene symbol.*/
    private HashMap<String,Gene> geneMap=null;
    /** Path to serialized Jannovar file with transcript information (e.g. ucsc.ser) */
    private String pathToSerializedJannovarFile=null;
    /** Path to a Variant Call Format (VCF) file that we will analyse. */
    private String pathToVCFfile=null;
    /** Path to a PED file (pedigree file) representing the samples on the VCF file. */
    private String pathToPedFile=null;
    /** Name of the outfile (default: jped.txt). */
    private String outfilename="jped.txt";
    /** An enumeration to indicate the mode inheritance, with AR=autosomal recessive, 
     * AD=autosomal dominant, X=X chromosomal, and HOM=autosomal recessive (only homozygous),
     * and COMPHET=autosomal recessive (only compound heterozygous). */
    private enum inheritanceMode {AR, AD, X, HOM, COMPHET};
    /** Mode of inheritance */
    private inheritanceMode mode;
   


    /**
     * An inner class desiged to hold all variants that correspond to a given gene.
     */
    public class Gene {
	public String symbol;
	public ArrayList<String> annots;
	public ArrayList<Variant> vars;
	public Gene (String s) {
	    this.symbol = s;
	    this.vars = new ArrayList<Variant>();
	    this.annots = new ArrayList<String>();
	}
	public void addVariant(Variant v, String a) {
	    this.vars.add(v);
	    this.annots.add(a);
	}
	public String getAnnotations() {
	    StringBuilder sb = new StringBuilder();
	    sb.append(symbol + "\n");
	  for (String a : annots) {
	      sb.append(a + "\n");
	  }
	  return sb.toString();
	}
    }

    public JPed(String argv[]) {
	parseCommandLineArguments(argv);
	this.geneMap = new HashMap<String,Gene>();
    }

    public static void main(String argv[]) {
	try {
	    JPed jped = new JPed(argv);
	    jped.deserializeUCSCdata();
	    jped.annotateVCF();
	    jped.processVariants();
	    jped.parsePedFile();
	    jped.filterByInheritance();
	    
	} catch (Exception e) {
	    e.printStackTrace();
	    JPed.usage();
	    System.exit(1);
	}
    }

    
    /**
     * This function is responsible for the input of the serialized
     * transcript data file that is produced by the main Jannovar
     * application and contains information about all transcripts
     * The transcript models are put into 
     * {@link #chromosomeMap}, which will allow efficient searching.
     */
    public void deserializeUCSCdata() throws JannovarException {
	SerializationManager manager = new SerializationManager();
	ArrayList<TranscriptModel> kgList = manager.deserializeKnownGeneList(this.pathToSerializedJannovarFile);
	this.chromosomeMap = Chromosome.constructChromosomeMapWithIntervalTree(kgList);
    }  
    
    /**
     * This function inputs a VCF file, and creates a list of all
     * variants identified ({@link #variantList}) as well as of
     * all samples names ({@link #sampleNames}).
     */
    public void annotateVCF() throws JannovarException {
	VCFReader parser = new VCFReader(this.pathToVCFfile);
	parser.parseFile();
	this.variantList = parser.getVariantList();
	this.sampleNames = parser.getSampleNames();
    }

   
    
    /**
     * For this example, we group the variants according to 
     * the gene they correspond to. In a real program, you 
     * would implement the logic of your analysis in a function(s)
     * like this, and filter variants according to your criteria
     * of interest. For instance, you might remove variants
     * that are common in the population by comparing them with 
     * data from dbSNP, or you might prioritize the variants by
     * checking their function with GO data or comparing the
     * phenotype of your patient with HPO data.
     */
    public void processVariants() throws JannovarException {
	int c=0;
	for (Variant v: this.variantList) {
	    String annot = annotateVariant(v);
	    if (annot==null)
		continue; /* Not found or not high enough priority */
	    if (v.isOffExomeTarget()) {
		continue; /* skip intronic and intergenic */
	    }
	    String sym = v.getGeneSymbol();
	    Gene g = null;
	    if (this.geneMap.containsKey(sym)) {
		  g = this.geneMap.get(sym);
	    } else {
		c++;
		if (c%500==0) {
		    System.out.println("Gene " + c+ ": " + sym);
		}
		 g = new Gene(sym);
		 this.geneMap.put(sym,g);
	     }
	     g.addVariant(v,annot);
	}
    }
   
      /**
     * This function outputs a single line in Jannovar format.
     * @param v The current variant with one or more annotations
     * @return A string with the HGVS variant nomenclature.
     */
    private String annotateVariant(Variant v) throws JannovarException
    {
	byte chr =  v.getChromosomeAsByte();
	String chrStr = v.get_chromosome_as_string();
	int pos = v.get_position();
	String ref = v.get_ref();
	String alt = v.get_alt();
	String gtype = v.getGenotypeAsString();
	Chromosome c = this.chromosomeMap.get(chr);
	AnnotationList annL = c.getAnnotationList(pos,ref,alt);
	v.setAnnotation(annL);
	String annotation = annL.getSingleTranscriptAnnotation();
	String effect = annL.getVariantType().toString();
	int prioritylevel = VariantType.priorityLevel(annL.getVariantType());
	if (prioritylevel > 1)
	    return null;
	String sym = annL.getGeneSymbol();
	String s = String.format("%s\t%s\t%s\t%s\t%s",
				 effect,sym,annotation,v.get_chromosomal_variant(),gtype);
	return s;
	
    }

    
    /**
     * Inputs a Ped file for multi-sample VCF files. If necessary,
     * also adjust the order of samples in the Pedigree to match the
     * order in the VCF file.
     */
    public void parsePedFile() throws JannovarException {
	PedFileParser parser = new PedFileParser();
	this.pedigree = parser.parseFile(this.pathToPedFile);
	this.pedigree.adjustSampleOrderInPedFile(this.sampleNames);
    }
    


    public void filterByInheritance() {
	try{
	    FileWriter fstream = new FileWriter(this.outfilename);
	    BufferedWriter out = new BufferedWriter(fstream);
	    out.write("The following genes have variants that are compatible with ");
	    
	    switch (this.mode) {
	    case AR:
		out.write("autosomal recessive inheritance\n");
		break;
	    case HOM:
		out.write("autosomal recessive inheritance (homozygous variants only)\n");
		break;
	    case COMPHET:
		out.write("autosomal recessive inheritance (compound heterozygous variants only)\n");
		break;
	    case AD:
		out.write("autosomal dominant inheritance\n");
		break;
	    case X:
		out.write("X chromosomal inheritance\n");
		break;
	    }
	    out.write(this.pedigree.getPedigreeSummary()+"\n");
	    for (Gene g : this.geneMap.values()) {
		ArrayList<Variant> varList = g.vars;
		if (this.mode==inheritanceMode.AR && 
		    this.pedigree.isCompatibleWithAutosomalRecessive(varList)) {
		    out.write(g.getAnnotations());
		} 
		else if (this.mode==inheritanceMode.HOM && 
		    this.pedigree.isCompatibleWithAutosomalRecessiveHomozygous(varList)) {
		    out.write(g.getAnnotations());
		} 
		else if (this.mode==inheritanceMode.COMPHET && 
		    this.pedigree.isCompatibleWithAutosomalRecessiveCompoundHet(varList)) {
		    out.write(g.getAnnotations());
		} 
		else if (this.mode==inheritanceMode.AD &&
		    this.pedigree.isCompatibleWithAutosomalDominant(varList)) {
		    out.write(g.getAnnotations());
		} 
		else if (this.mode==inheritanceMode.X &&
		    this.pedigree.isCompatibleWithXChromosomalRecessive(varList)) {
		     out.write(g.getAnnotations());
		} 
	    }
	    out.close();
	} catch (IOException e) {
	    e.printStackTrace();
	}
    } 

 

     /**
     * Parse the command line using apache's CLI. A copy of the library is included
     * in the Jannovar tutorial archive in the lib directory.
     * @see "http://commons.apache.org/proper/commons-cli/index.html"
     */
    private void parseCommandLineArguments(String[] args)
    {
	try {
	    Options options = new Options();
	    options.addOption(new Option("H","help",false,"Shows this help"));
	    options.addOption(new Option("V","vcf",true,"Path to VCF file."));
	    options.addOption(new Option("D","deserialize",true,"Path to jannovar transcript definition file (e.g., ucsc.ser)."));
	    options.addOption(new Option("P","ped",true,
					 "Path to ped file"));
	    options.addOption(new Option("I","inheritance",true,
					 "Inheritance pattern (AR, AD, X)"));
	    options.addOption(new Option("O","outfile",true,
					  "name of outfile (default: jlink.txt)"));
	    Parser parser = new GnuParser();
	    CommandLine cmd = parser.parse(options,args);
	    if ( cmd.hasOption("H")){
		usage();
	    }
	    if (cmd.hasOption("D")) {
		this.pathToSerializedJannovarFile = cmd.getOptionValue("D");
	    } else {
		JPed.usage();
	    }
	    if (cmd.hasOption("V")) {
		this.pathToVCFfile = cmd.getOptionValue("V");
	    } else {
		usage();
	    }
	    if (cmd.hasOption("P")) {
		this.pathToPedFile = cmd.getOptionValue("P");
	    } else {
		usage();
	    }
	    if (cmd.hasOption("I")) {
		setModeOfInheritance(cmd.getOptionValue("I"));
	    } else {
		usage();
	    }
	    if (cmd.hasOption("O")) {
		this.outfilename = cmd.getOptionValue("O");
	    } 
	} catch (ParseException pe) {
	    System.err.println("Error parsing command line options");
	    System.err.println(pe.getMessage());
	    System.exit(1);
	}
    }

    
    /**
     * Parse a command line argument representing
     * the inheritance mode (should be one of
     * AR, AD, or X; otherwise the program will
     * terminate with an error).
     */
    private void setModeOfInheritance(String moi) {
	if (moi.equalsIgnoreCase("AR"))
	    this.mode = inheritanceMode.AR;
	else if (moi.equalsIgnoreCase("AD"))
	    this.mode = inheritanceMode.AD;
	else if (moi.equalsIgnoreCase("X"))
	    this.mode = inheritanceMode.X;
	else if (moi.equalsIgnoreCase("AR-HOM"))
	    this.mode = inheritanceMode.HOM;
	else if (moi.equalsIgnoreCase("AR-COMPHET"))
	    this.mode = inheritanceMode.COMPHET;
	else {
	    System.out.println("Did not recognize mode of inheritance: \"" + moi + "\"");
	    System.out.println("Please use one of \"AR\", \"AR-HOM\", \"AR-COMPHET\", \"AD\", or \"X\", e.g., \"-I AD\"");
	    System.exit(1);
	}
    }


    /**
     * Prints a usage message to teh console.
     */
    public static void usage() {
	System.err.println("[INFO] Usage: java -jar JPed.jar -D ucsc_hg19.ser -V sample.vcf -P sample.ped  -I AR [-O fname]");
	System.err.println("[INFO] where");
	System.err.println("[INFO]");
	System.err.println("[INFO] -D: the serialized transcript data file from Jannovar (e.g., ucsc_hg19.ser)");
	System.err.println("[INFO] -V: the VCF file representing samples from a family");
	System.err.println("[INFO] -P: the corresponding PED file");
	System.err.println("[INFO] -I: the mode of inheritance to filter for (AD, AR, X, AR-HOM, AR-COMPHET)");
    	System.err.println("[INFO] -O: Name of output file (optional)");
	System.exit(1);
    }
   
}
