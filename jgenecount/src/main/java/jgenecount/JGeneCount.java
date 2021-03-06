package jgenecount;

import java.util.ArrayList;
import java.util.HashSet;

/** Command line parser from apache */
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.GnuParser;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.commons.cli.Parser;

import bsh.org.objectweb.asm.Constants;

import com.google.common.collect.ImmutableMap;

import de.charite.compbio.jannovar.JannovarException;
import de.charite.compbio.jannovar.annotation.VariantAnnotator;
import de.charite.compbio.jannovar.data.Chromosome;
import de.charite.compbio.jannovar.data.JannovarData;
import de.charite.compbio.jannovar.data.JannovarDataSerializer;
import de.charite.compbio.jannovar.data.ReferenceDictionary;
import de.charite.compbio.jannovar.htsjdk.VariantContextAnnotator;
import de.charite.compbio.jannovar.reference.TranscriptModel;
/** Classes from the Jannovar library. */

/**
 * This application counts and lists all genes located in the chromosomal interval that is passed via the command line.
 * For instance, the command <code>$ java -jar JGeneCount.jar -D ucsc_hg19.ser -I chr14:29781404-30552936</code> will
 * yield:
 *
 * <pre>
 * Genes and transcripts in 14:29781404-30552936
 * 1) MIR548AI
 * 2) PRKD1
 * 3) BC062469
 * 4) U6
 * 4 genes with 1 coding and 4 noncoding transcripts
 * </pre>
 *
 * This kind of information may be useful when evaluating copy number variants, for instance.
 *
 * @author Peter Robinson
 * @version 0.1 (January 19,2014)
 */
public class JGeneCount implements Constants {

	/** Path to serialized Jannovar file with transcript information (e.g. ucsc.ser) */
	private String pathToSerializedJannovarFile = null;
	/** Name of the outfile (default: jped.txt). */
	private String outfilename = "jgenecount.txt";
	/** Map of Chromosomes (Jannovar object storing all transcripts from the Chromosome). */
	private ImmutableMap<Integer, Chromosome> chromosomeMap = null;
	/** Current chromosome ID. */
	private Integer chromosome = null;
	/** Genome chromosome information */
	private ReferenceDictionary refDict = null;
	/** {@link VariantAnnotator} to use for generating variant annotations */
	private VariantContextAnnotator annotator = null;
	/** Start position of the interval on the chromosome */
	private int from;
	/** End position of the interval on the chromosome */
	private int to;
	/** List of all transcripts found in the interval */
	ArrayList<TranscriptModel> transcriptsInInterval = null;

	public JGeneCount(String args[]) {
		parseCommandLineArguments(args);
	}

	public static void main(String[] args) {
		JGeneCount jgc = new JGeneCount(args);
		try {
			jgc.deserializeUCSCdata();
		} catch (JannovarException e) {
			e.printStackTrace();
			System.exit(1);
		}
		jgc.printGenesInInterval();
	}

	/**
	 * Prints a list of the genes in the interval to standard out.
	 */
	public void printGenesInInterval() {
		HashSet<String> symbolSet = new HashSet<String>();
		int n_coding = 0;
		int n_noncoding = 0;
		for (TranscriptModel tm : this.transcriptsInInterval) {
			symbolSet.add(tm.getGeneSymbol());
			if (tm.isCoding())
				n_coding++;
			else
				n_noncoding++;
		}
		System.out.println("Genes and transcripts in " + this.chromosome + ":" + this.from + "-" + this.to);
		int c = 0;
		for (String s : symbolSet) {
			c++;
			System.out.println(c + ") " + s);
		}
		System.out.println(symbolSet.size() + " genes with " + n_coding + " coding and " + n_noncoding
				+ " noncoding transcripts");
	}

	/**
	 * Parses an input string such as chr4:123-456 into the chromosome, from, and to parts.
	 */
	private void setInterval(String interval) {
		String A[] = interval.split(":");
		if (A.length != 2) {
			System.err.println("Malformed interval string: " + interval);
			System.err.println("Require string such as \"chr1:123-456\"");
			System.exit(1);
		}

		this.chromosome = refDict.getContigNameToID().get(A[0]).intValue();

		String B[] = A[1].split("-");
		if (B.length != 2) {
			System.err.println("Malformed interval string: " + interval);
			System.err.println("Require string such as \"chr1:123-456\"");
			System.exit(1);
		}
		try {
			this.from = Integer.parseInt(B[0]);
			this.to = Integer.parseInt(B[1]);
		} catch (NumberFormatException E) {
			System.err.println("Malformed intervalstring: " + A[1]);
			System.err.println("Require string such as \"chr1:123-456\"");
			System.err.print(E);
			System.exit(1);
		}
	}

	/**
	 * This function is responsible for the input of the serialized transcript data file that is produced by the main
	 * Jannovar application. It places all transcripts that are located in the chromosomal interval into
	 * {@link #transcriptsInInterval}, and discards the remaining ones.
	 */
	public void deserializeUCSCdata() throws JannovarException {
		JannovarData data = new JannovarDataSerializer(this.pathToSerializedJannovarFile).load();
		this.chromosomeMap = data.getChromosomes();
		this.refDict = data.getRefDict();
		this.annotator = new VariantContextAnnotator(data.getRefDict(), data.getChromosomes());
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

			options.addOption(new Option("D", "deserialize", true,
					"Path to jannovar transcript definition file (e.g., ucsc.ser)."));
			options.addOption(new Option("I", "interval", true, "Chromosomeal interval, e.g., chr2:23456-8172634"));
			Parser parser = new GnuParser();
			CommandLine cmd = parser.parse(options, args);
			if (cmd.hasOption("H")) {
				usage();
			}
			if (cmd.hasOption("D")) {
				this.pathToSerializedJannovarFile = cmd.getOptionValue("D");
			} else {
				JGeneCount.usage();
			}
			if (cmd.hasOption("I")) {
				setInterval(cmd.getOptionValue("I"));
			} else {
				JGeneCount.usage();
			}
		} catch (ParseException pe) {
			System.err.println("Error parsing command line options");
			System.err.println(pe.getMessage());
			System.exit(1);
		}
	}

	/**
	 * Prints a usage message to teh console.
	 */
	public static void usage() {
		System.err.println("[INFO] Usage: java -jar JGeneCount.jar -D ucsc_hg19.ser -I chr4:123-456");
		System.err.println("[INFO] where");
		System.err.println("[INFO]");
		System.err.println("[INFO] -D: the serialized transcript data file from Jannovar (e.g., ucsc_hg19.ser)");
		System.err
				.println("[INFO] -I: the chromosomal interval  (for instance chr4:1223-456, from position 123 to 456 on chromosome 4)");
		System.exit(1);
	}

}