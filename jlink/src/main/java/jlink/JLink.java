package jlink;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.GnuParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.commons.cli.Parser;

import com.google.common.base.Joiner;
import com.google.common.collect.ImmutableMap;

import de.charite.compbio.jannovar.JannovarException;
import de.charite.compbio.jannovar.annotation.Annotation;
import de.charite.compbio.jannovar.annotation.PutativeImpact;
import de.charite.compbio.jannovar.annotation.VariantAnnotations;
import de.charite.compbio.jannovar.annotation.VariantAnnotator;
import de.charite.compbio.jannovar.data.Chromosome;
import de.charite.compbio.jannovar.data.JannovarData;
import de.charite.compbio.jannovar.data.JannovarDataSerializer;
import de.charite.compbio.jannovar.data.ReferenceDictionary;
import de.charite.compbio.jannovar.htsjdk.InvalidCoordinatesException;
import de.charite.compbio.jannovar.htsjdk.VariantContextAnnotator;

/** Command line parser from apache */

/**
 * Example Jannovar program designed to demonstrate how to perform filtering of a VCF file acccording to linkage
 * interval. filtering with exome variant annotation.
 *
 * @author Peter Robinson
 * @version 1.0 (December 17, 2013)
 */
public class JLink {
	/** Map of Chromosomes (Jannovar object storing all transcripts from the Chromosome). */
	private ImmutableMap<Integer, Chromosome> chromosomeMap = null;
	/** Genome chromosome information */
	private ReferenceDictionary refDict = null;
	/** {@link VariantAnnotator} to use for generating variant annotations */
	private VariantContextAnnotator annotator = null;
	/** List of variants from input file to be analysed. */
	private ArrayList<VariantContext> variantList = null;
	/** The chromosome of the linkage interval we are filtering on. */
	private int chromosome;
	/** The start position of the linkage interval. */
	private int fromPosition;
	/** The start position of the linkage interval. */
	private int toPosition;
	/** Name of the outfile (default: jlink.txt). */
	private String outfilename = "jlink.txt";

	/** Path to serialized Jannovar file with transcript information (e.g. ucsc.ser) */
	private String pathToSerializedJannovarFile = null;
	/** Path to a Variant Call Format (VCF) file that we will analyse. */
	private String pathToVCFfile = null;
	/** String of the form "chr4:12345-54321" representing a libnkage interval */
	private String intervalString = null;

	public static void main(String argv[]) throws InvalidCoordinatesException {
		JLink jlink = new JLink(argv);
		try {
			jlink.deserializeUCSCdata();
			jlink.annotateVCF();
			jlink.filterVariants();
		} catch (JannovarException e) {
			e.printStackTrace();
		} catch (IOException ioe) {
			System.out.println("Could not open file for output");
			ioe.printStackTrace();
		}
	}

	/**
	 * The constructor reads the command line arguments and initializes the linkage interval to filter over.
	 */
	public JLink(String argv[]) {
		parseCommandLineArguments(argv);
		parseInterval();
	}

	/**
	 * This function filters out variants that are not in the linkage interval and prints out the variants that are in
	 * the interval to file.
	 *
	 * @throws InvalidCoordinatesException
	 *             on problems with coordiante computations
	 */
	public void filterVariants() throws JannovarException, IOException, InvalidCoordinatesException {
		FileWriter fstream = new FileWriter(this.outfilename);
		BufferedWriter out = new BufferedWriter(fstream);
		for (VariantContext vc : this.variantList) {
			int chrID = refDict.getContigNameToID().get(vc.getChr());
			if (chrID != this.chromosome)
				continue;
			if (vc.getStart() < this.fromPosition)
				continue;
			if (vc.getEnd() > this.toPosition)
				continue;
			String a = annotateVariant(vc);
			if (a == null)
				continue;
			out.write(a + "\n");
		}
		out.close();
	}

	/**
	 * This function outputs a single line in Jannovar format.
	 *
	 * @param vc
	 *            The current variant with one or more annotations
	 * @return A string with the HGVS variant nomenclature.
	 * @throws InvalidCoordinatesException
	 *             on problems with coordiante computations
	 */
	private String annotateVariant(VariantContext vc) throws JannovarException, InvalidCoordinatesException {
		// TODO(holtgrewe): Nick, we now have multiple alternative alleles and annotations, please check logic.
		// TODO(holtgrewe): Also, we only look at the first genotype here.
		for (VariantAnnotations annos : annotator.buildAnnotations(vc)) {
			for (Annotation anno : annos.getAnnotations()) {
				if (anno.getPutativeImpact() != PutativeImpact.MODERATE)
					continue; // skip low-level impacts
				return Joiner.on('\t').join(anno.getEffects().first(), anno.getGeneSymbol(), anno.getGenomeVariant(),
						vc.getGenotype(0));
			}
		}
		return null; // no valid annotation
	}

	/**
	 * Parse a string such as chr4:12345-54321 into the chromosome, start, and stop positions, and use them to
	 * initialize the class variables {@link #chromosome}, {@link #fromPosition}, and {@link #toPosition}.
	 */
	public void parseInterval() {
		int idx = this.intervalString.indexOf(":");
		if (idx < 0) {
			System.out.println("Malformed interval string: " + this.intervalString);
			System.exit(1);
		}
		String chr = this.intervalString.substring(0, idx);
		String rest = this.intervalString.substring(idx + 1);
		idx = rest.indexOf("-");
		if (idx < 0) {
			System.out.println("Malformed interval string: " + this.intervalString);
			System.exit(1);
		}
		String from = rest.substring(0, idx);
		String to = rest.substring(idx + 1);
		this.chromosome = refDict.getContigNameToID().get(chr).intValue();
		this.fromPosition = Integer.parseInt(from);
		this.toPosition = Integer.parseInt(to);
	}

	/**
	 * Users of JLink must supply the path to a serialized transcript information file as created by Jannovar. This
	 * function deserializes this file.
	 */
	public void deserializeUCSCdata() throws JannovarException {
		JannovarData data = new JannovarDataSerializer(this.pathToSerializedJannovarFile).load();
		this.chromosomeMap = data.getChromosomes();
		this.refDict = data.getRefDict();
		this.annotator = new VariantContextAnnotator(data.getRefDict(), data.getChromosomes());
	}

	/**
	 * This function inputs a VCF file, and prints the annotated version thereof to a file (name of the original file
	 * with the suffix .jannovar).
	 */
	public void annotateVCF() throws JannovarException {
		VCFFileReader parser = new VCFFileReader(new File(this.pathToVCFfile), /* requireIndex= */false);
		this.variantList = new ArrayList<VariantContext>();
		for (VariantContext vc : parser)
			this.variantList.add(vc);
		parser.close();
	}

	/**
	 * Parse the command line using apache's CLI. A copy of the library is included in the Jannovar tutorial archive in
	 * the lib directory.
	 *
	 * @see "http://commons.apache.org/proper/commons-cli/index.html"
	 */
	private void parseCommandLineArguments(String[] args) {
		try {
			Options options = new Options();
			options.addOption(new Option("H", "help", false, "Shows this help"));
			options.addOption(new Option("V", "vcf", true, "Path to VCF file."));
			options.addOption(new Option("D", "deserialize", true,
					"Path to jannovar transcript definition file (e.g., ucsc.ser)."));
			options.addOption(new Option("I", "interval", true, "interval to filter on (e.g., chr4:12345-54321)"));
			options.addOption(new Option("O", "outfile", true, "name of outfile (default: jlink.txt)"));
			Parser parser = new GnuParser();
			CommandLine cmd = parser.parse(options, args);
			if (cmd.hasOption("H")) {
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
