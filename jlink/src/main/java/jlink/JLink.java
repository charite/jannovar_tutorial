package jlink;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;
import java.util.Collections;

import java.io.FileWriter; 
import java.io.BufferedWriter;
import java.io.IOException;

/** Command line parser from apache */
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.GnuParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.commons.cli.Parser;

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
 * Example Jannovar program designed to demonstrate how to perform
 * filtering of a VCF file acccording to linkage interval.
 * filtering with exome variant annotation.
 * @author Peter Robinson
 * @version 1.0 (December 17, 2013)
 */
public class JLink {
    /** Map of Chromosomes */
    private HashMap<Byte,Chromosome> chromosomeMap=null;
     /** List of variants from input file to be analysed. */
    private ArrayList<Variant> variantList=null;
    /** List of all sample names of VCF file */
    static int X_CHROMOSOME =23;
    static int Y_CHROMOSOME =24;
    static int M_CHROMOSOME =25;
    /** The chromosome of the linkage interval we are filtering on. */
    private int chromosome;
    /** The start position of the linkage interval. */
    private int fromPosition;
    /** The start position of the linkage interval. */
    private int toPosition;
    /** Name of the outfile (default: jlink.txt). */
    private String outfilename="jlink.txt";

    /** Path to serialized Jannovar file with transcript information (e.g. ucsc.ser) */
    private String pathToSerializedJannovarFile=null;
    /** Path to a Variant Call Format (VCF) file that we will analyse. */
    private String pathToVCFfile=null;
    /** String of the form "chr4:12345-54321" representing a libnkage interval */
    private String intervalString=null;

    public static void main(String argv[]) {
	  JLink jlink = new JLink(argv);
	  try {
	      jlink.deserializeUCSCdata();
	      jlink.annotateVCF();
	      jlink.filterVariants();
	  } catch (JannovarException e) {
	    e.printStackTrace();
	  } catch(IOException ioe) {
	      System.out.println("Could not open file for output");
	      ioe.printStackTrace();
	  }    
    }

    /**
     * The constructor reads the command line arguments and initializes the
     * linkage interval to filter over.
     */
    public JLink(String argv[]) {
	parseCommandLineArguments(argv);
	parseInterval();
    }


    /**
     * This function filters out variants that are not in the linkage interval
     * and prints out the variants that are in the interval to file.
     */
    public void filterVariants() throws JannovarException, IOException  {
	FileWriter fstream = new FileWriter(this.outfilename);
	BufferedWriter out = new BufferedWriter(fstream);
	for (Variant v: this.variantList) {
	    if (v.get_chromosome()  != this.chromosome)
		continue;
	    if (v.get_position() < this.fromPosition)
		continue;
	    if (v.get_position() > this.toPosition)
		continue;
	    String a = annotateVariant(v);
	    if (a == null)
		continue;
	    out.write(a + "\n");
	}
	out.close();
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
	VariantType vt = annL.getVariantType();
	int level = VariantType.priorityLevel(vt);
	 if (level > 1)
	     return null;
	String effect = annL.getVariantType().toString();
	String sym = annL.getGeneSymbol();
	String s = String.format("%s\t%s\t%s\t%s\t%s",
				 effect,sym,annotation,v.get_chromosomal_variant(),gtype);
	return s;
	
    }


    /**
     * Parse a string such as chr4:12345-54321 into the chromosome, start, and stop
     * positions, and use them to initialize the class variables
     * {@link #chromosome}, {@link #fromPosition}, and {@link #toPosition}.
     */
    public void parseInterval() {
	int idx = this.intervalString.indexOf(":");
	if (idx < 0) {
	    System.out.println("Malformed interval string: " + this.intervalString);
	    System.exit(1);
	}
	String chr = this.intervalString.substring(0,idx);
	String rest = this.intervalString.substring(idx+1);
	idx = rest.indexOf("-");
	if (idx < 0) {
	    System.out.println("Malformed interval string: " + this.intervalString);
	    System.exit(1);
	}
	String from = rest.substring(0,idx);
	String to = rest.substring(idx+1);
	if (chr.startsWith("chr"))
	    chr = chr.substring(3);
	else {
	    System.out.println("Malformed interval string: " + this.intervalString);
	    System.exit(1);
	}
	if (chr.equals("X")) {
	    this.chromosome = X_CHROMOSOME ;
	} else if (chr.equals("Y")) {
	 this.chromosome = Y_CHROMOSOME ;
	}else if (chr.equals("M")) {
	 this.chromosome = M_CHROMOSOME ;
	} else {
	    this.chromosome = Integer.parseInt(chr);
	}
	this.fromPosition = Integer.parseInt(from);
	this.toPosition = Integer.parseInt(to);
    }


    /**
     * Users of JLink must supply the path to a serialized transcript information
     * file as created by Jannovar. This function deserializes this file. */      
    public void deserializeUCSCdata() throws JannovarException  {
	SerializationManager manager = new SerializationManager();
	ArrayList<TranscriptModel> kgList = manager.deserializeKnownGeneList(this.pathToSerializedJannovarFile);
	this.chromosomeMap = Chromosome.constructChromosomeMapWithIntervalTree(kgList);

    }  
    
    /**
     * This function inputs a VCF file, and prints the annotated version thereof
     * to a file (name of the original file with the suffix .jannovar).
     */
    public void annotateVCF() throws JannovarException {
	VCFReader parser = new VCFReader(this.pathToVCFfile);
	parser.parseFile();
	this.variantList = parser.getVariantList();
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
	    options.addOption(new Option("I","interval",true,
					 "interval to filter on (e.g., chr4:12345-54321)"));
	    options.addOption(new Option("O","outfile",true,
					  "name of outfile (default: jlink.txt)"));
	    Parser parser = new GnuParser();
	    CommandLine cmd = parser.parse(options,args);
	    if ( cmd.hasOption("H")){
		usage(options);
	    }
	    if (cmd.hasOption("D")) {
		this.pathToSerializedJannovarFile = cmd.getOptionValue("D");
	    } else {
		usage(options);
	    }
	    if (cmd.hasOption("V")) {
		this.pathToVCFfile = cmd.getOptionValue("V");
	    } else {
		usage(options);
	    }
	    if (cmd.hasOption("I")) {
		this.intervalString = cmd.getOptionValue("I");
	    } else {
		usage(options);
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

    private void usage(Options options) {
	HelpFormatter formatter = new HelpFormatter();
	String helpString = "java -jar JLink.jar -V vcffile -D jfile -I interval [-O ofile]";
	formatter.printHelp(helpString, options);
	System.exit(0);
    }

    



}
