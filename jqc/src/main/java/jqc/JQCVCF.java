package jqc;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.GnuParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.commons.cli.Parser;

import de.charite.compbio.jannovar.JannovarException;
import de.charite.compbio.jannovar.data.Chromosome;
import de.charite.compbio.jannovar.reference.TranscriptModel;
/** Command line parser from apache */

/**
 * Example Jannovar program designed to demonstrate how to extract some Q/C statistics for an exome or other VCF file.
 *
 * @author Peter Robinson
 * @version 1.2 (December 29, 2013)
 */
public class JQCVCF {
	/** Map of Chromosomes */
	private HashMap<Byte, Chromosome> chromosomeMap = null;
	/** List of variants from input file to be analysed. */
	private ArrayList<Variant> variantList = null;
	/** The names of the samples in the VCF file (they are parsed automatically from the file). */
	ArrayList<String> sampleNames = null;
	/** Name of the outfile (default: jqc.txt). */
	private String outfilename = "jqc.tex";
	/** Name of the original VCF file (to be set from the command line). */
	private String vcfFileName = null;
	/** Path to serilaized Jannovar file with transcript information (e.g. ucsc_hg19.ser) */
	private String pathToSerializedJannovarFile = null;
	/** Path to a Variant Call Format (VCF) file that we will analyse. */
	private String pathToVCFfile = null;
	/** Total number of called variants . */
	private int n_totalVariants = 0;
	/**
	 * Flag whether the program should attempt to automatically compile a PDF file as output. This will work on most
	 * linux systems but probably not on Windows or Mac (You can adjust the code to match the tools and paths on your
	 * system if desired. See {@link jqc.LatexReportWriter LatexReportWriter}.
	 */
	private boolean compilePDF = false;
	/** Number of samples in the VCF file. */
	private int n_samples = 0;
	/** A statistic object to calculate stats for the PHRED variant quality score */
	private DescriptiveStatistic phredStat = null;
	/**
	 * Object that will write the latex file with the Q/C report and that can attempt to compile this file to PDF.
	 */
	private LatexReportWriter writer = null;
	/** One {@link SampleReport} object is created for each sample in the VCF file. */
	private ArrayList<SampleReport> samples = null;
	/**
	 * A flag to indicate whether we should disregard variants that are more than {@link #distanceThreshold} removed
	 * from a target exon.
	 */
	private boolean doTargetFilter = false;
	/**
	 * A distance threshold. If a variant is more nucleotides distance from the target region (any exon), it will not be
	 * included in the calculation of quality quantities such as tv/ti. Instead, we will just record the number of
	 * offTarget variants.
	 */
	private int distanceThreshold;
	/** Number of off target variants (only calculated if {@link #doTargetFilter} is true). */
	private int n_offTarget = 0;
	/**
	 * Upper maximum limit for the Phred score of 10,000. Occasioanlly, a value of such as 6442456297.26 gets reported,
	 * and we want to cap that.
	 */
	private float THRESHOLD_PHRED = 10000f;

	public static void main(String argv[]) {
		JQCVCF jqc = new JQCVCF(argv);
		try {
			jqc.deserializeUCSCdata();
			jqc.inputVCF();
			jqc.outputLatexReport();
			if (jqc.doCompilePDF()) {
				jqc.compilePDF();
			}
		} catch (JannovarException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	/**
	 * The constructor initializes some variables based on the command line arguments.
	 *
	 * @param argv
	 *            Command line arguments.
	 */
	public JQCVCF(String argv[]) {
		parseCommandLineArguments(argv);
		this.samples = new ArrayList<SampleReport>();
		this.phredStat = new DescriptiveStatistic("PHRED Score");
	}

	/**
	 * This function inputs a VCF file, and prints the annotated version thereof to a file (name of the original file
	 * with the suffix .jannovar).
	 */
	public void inputVCF() throws JannovarException {
		VCFReader parser = new VCFReader(this.pathToVCFfile);
		parser.inputVCFheader();
		this.n_samples = parser.getNumberOfSamples();
		this.sampleNames = parser.getSampleNames();
		/**
		 * Initiliaze the SampleReport objects, one for each sample in the VCF file.
		 */
		int index = 0;
		for (String name : this.sampleNames) {
			SampleReport report = new SampleReport(name, index);
			index++;
			this.samples.add(report);
		}
		this.vcfFileName = parser.getVCFFileName();
		/* Iterate over all variants in the VCF file. */
		Iterator<Variant> it = parser.getVariantIterator();

		while (it.hasNext()) {
			Variant v = it.next();
			byte chr = v.getChromosomeAsByte();
			Chromosome c = this.chromosomeMap.get(chr);
			v.annotate(c);
			this.n_totalVariants++;
			boolean ontarget = doQualityControl(v);
			if (ontarget) {
				float phred = v.getVariantPhredScore();
				phred = Math.min(phred, THRESHOLD_PHRED);
				this.phredStat.addValue(phred);
			}
		}
		// Now add the phred results to each of the reports
		this.phredStat.evaluate();
		for (SampleReport r : this.samples) {
			r.setPhredString(this.phredStat.toString());
		}
	}

	private boolean isOffTarget(Variant v) {
		VariantType vt = v.getVariantTypeConstant();
		String ann = v.getAnnotation();
		if (vt == VariantType.INTERGENIC)
			return true;
		if (vt == VariantType.UPSTREAM || vt == VariantType.DOWNSTREAM) {
			// e.g., KRTAP19-2(dist=10)
			int i = ann.indexOf("=");
			int j = ann.indexOf(")");
			if (i < 0 || j < 0) {
				System.err.println("[WARNING] Malformed Upstream annotation");
				System.err.println("[WARNING] " + ann);
				return true;
			}
			Integer ii = Integer.parseInt(ann.substring(i + 1, j));
			if (ii > this.distanceThreshold)
				return true;
			else
				return false;
		}
		if (vt == VariantType.INTRONIC || vt == VariantType.ncRNA_INTRONIC) {
			// e.g. TMBIM6(uc001rux.2:dist to exon4=2439;dist to exon5=33)
			int i = ann.indexOf("=");
			int j = ann.indexOf(";");
			Integer ii = 0;
			try {
				ii = Integer.parseInt(ann.substring(i + 1, j));
			} catch (NumberFormatException e) {
				System.err.println("[ERROR] Malformed annotation:" + ann);
				System.err.println("[ERROR] skipping this variant");
				return true;
			}
			if (ii <= this.distanceThreshold)
				return false;
			i = ann.indexOf("=", i + 1);
			j = ann.indexOf(")", i);
			try {
				ii = Integer.parseInt(ann.substring(i + 1, j));
			} catch (NumberFormatException e) {
				System.err.println("[ERROR] Malformed annotation:" + ann);
				System.err.println("[ERROR] skipping this variant");
				return true;
			}
			if (ii > this.distanceThreshold)
				return true;
			else
				return false;
		}
		return false;
	}

	/**
	 * This function adds the Variant and its quality parameters to each of the SampleReports (one for each sample in
	 * the VCF file).
	 *
	 * @param v
	 *            The current Variant.
	 */
	private boolean doQualityControl(Variant v) {
		if (this.doTargetFilter && isOffTarget(v)) {
			this.n_offTarget++;
			return false;
		}
		for (int i = 0; i < this.samples.size(); ++i) {
			SampleReport rep = this.samples.get(i);
			rep.addVariant(v);
		}
		return true;
	}

	/**
	 * This function will use the runtime environment to call pdflatex from the shell in order to compile the latex file
	 * that was written by {@link #outputLatexReport}.
	 */
	public void compilePDF() {
		this.writer.createPDFFile(this.outfilename);
	}

	/**
	 * This will write the results to a latex file and try to compile the results to a file called jqc.pdf.
	 */
	public void outputLatexReport() throws IOException {
		this.writer = new LatexReportWriter(this.vcfFileName);
		FileWriter fstream = new FileWriter(this.outfilename);
		BufferedWriter out = new BufferedWriter(fstream);
		try {
			writer.writeHeader(out);
			if (this.doTargetFilter) {
				writer.outputOffTargetSummary(out, this.distanceThreshold, this.n_offTarget, this.n_totalVariants);
			}
			Iterator<SampleReport> it = this.samples.iterator();
			while (it.hasNext()) {
				SampleReport report = it.next();
				report.writeLatexSummary(out);
			}
			writer.writeFooter(out);
			out.close();
		} catch (IOException e) {
			System.err.println("[ERROR] Could not write latex file");
			e.printStackTrace();
			System.exit(1);
		}
	}

	public boolean doCompilePDF() {
		return this.compilePDF;
	}

	/**
	 * Users of JLink must supply the path to a serialized transcript information file as created by Jannovar. This
	 * function deserializes this file.
	 */
	public void deserializeUCSCdata() throws JannovarException {
		SerializationManager manager = new SerializationManager();
		ArrayList<TranscriptModel> kgList = manager.deserializeKnownGeneList(this.pathToSerializedJannovarFile);
		this.chromosomeMap = Chromosome.constructChromosomeMapWithIntervalTree(kgList);
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
			options.addOption(new Option(null, "pdf", false, "compile PDF file with pdflatex"));
			options.addOption(new Option("O", "outfile", true, "name of outfile (default: jqc.txt)"));
			options.addOption(new Option(null, "threshold", true, "threshold distance to target for Q/C analysis"));
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
			if (cmd.hasOption("O")) {
				this.outfilename = cmd.getOptionValue("O");
			}
			if (cmd.hasOption("pdf")) {
				this.compilePDF = true;
			} else {
				this.compilePDF = false;
			}
			if (cmd.hasOption("threshold")) {
				String t = cmd.getOptionValue("threshold");
				this.distanceThreshold = Integer.parseInt(t);
				this.doTargetFilter = true;
			} else {
				this.doTargetFilter = false;
			}
		} catch (ParseException pe) {
			System.err.println("Error parsing command line options");
			System.err.println(pe.getMessage());
			System.exit(1);
		}
	}

	private void usage(Options options) {
		HelpFormatter formatter = new HelpFormatter();
		String helpString = "java -jar JQCVCF.jar -V vcffile -D jfile  [-O ofile]";
		formatter.printHelp(helpString, options);
		System.exit(0);
	}

}