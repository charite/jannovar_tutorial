package jqc;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Iterator;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.GnuParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.commons.cli.Parser;

import com.google.common.collect.ImmutableMap;
import com.google.common.collect.ImmutableSortedSet;

import de.charite.compbio.jannovar.JannovarException;
import de.charite.compbio.jannovar.annotation.Annotation;
import de.charite.compbio.jannovar.annotation.VariantAnnotations;
import de.charite.compbio.jannovar.annotation.VariantEffect;
import de.charite.compbio.jannovar.data.Chromosome;
import de.charite.compbio.jannovar.data.JannovarData;
import de.charite.compbio.jannovar.data.JannovarDataSerializer;
import de.charite.compbio.jannovar.data.ReferenceDictionary;
import de.charite.compbio.jannovar.htsjdk.InvalidCoordinatesException;
/** Command line parser from apache */
import de.charite.compbio.jannovar.htsjdk.VariantContextAnnotator;

/**
 * Example Jannovar program designed to demonstrate how to extract some Q/C statistics for an exome or other VCF file.
 *
 * @author Peter Robinson
 * @version 1.2 (December 29, 2013)
 */
public class JQCVCF {
	/** Map of Chromosomes */
	private ImmutableMap<Integer, Chromosome> chromosomeMap = null;
	/** List of variants from input file to be analysed. */
	private ImmutableMap<Integer, Chromosome> variantList = null;
	/** Information about genome (chromosome lengths and names). */
	private ReferenceDictionary refDict = null;
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
	private double THRESHOLD_PHRED = 10000f;

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
		} catch (InvalidCoordinatesException e) {
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
	 *
	 * @throws InvalidCoordinatesException
	 *             in the case of problems in the variant annotation
	 */
	public void inputVCF() throws JannovarException, InvalidCoordinatesException {
		// open file using HTSJDK's VCFFileReader and get info from headers
		VCFFileReader parser = new VCFFileReader(new File(this.pathToVCFfile), /* requireIndex= */false);
		this.n_samples = parser.getFileHeader().getNGenotypeSamples();
		this.sampleNames = parser.getFileHeader().getSampleNamesInOrder();

		/**
		 * Initiliaze the SampleReport objects, one for each sample in the VCF file.
		 */
		int index = 0;
		for (String name : this.sampleNames) {
			SampleReport report = new SampleReport(name, index);
			index++;
			this.samples.add(report);
		}
		this.vcfFileName = this.pathToVCFfile;
		/* Iterate over all variants in the VCF file and build annotations */
		VariantContextAnnotator annotator = new VariantContextAnnotator(refDict, chromosomeMap);
		for (VariantContext vc : parser) {
			this.n_totalVariants += 1;

			// TODO(holtgrew): Nick, here the behaviour changed from considering the first to all alleles.
			for (VariantAnnotations annos : annotator.buildAnnotations(vc)) { // consider each alternative allele
				boolean onTarget = doQualityControl(vc, annos); // QC of annotations for one allele/variant
				if (onTarget) {
					double phred = vc.getPhredScaledQual();
					phred = Math.min(phred, THRESHOLD_PHRED);
					this.phredStat.addValue(phred);
				}
			}
		}

		// Now add the phred results to each of the reports
		this.phredStat.evaluate();
		for (SampleReport r : this.samples) {
			r.setPhredString(this.phredStat.toString());
		}
	}

	/**
	 * Perform QC of one annotated variant.
	 *
	 * There might be more than one annotation (one for each affected transcript), thus the parameter is a
	 * {@link VariantAnnotations}.
	 */
	private boolean isOffTarget(VariantAnnotations va) {
		Annotation anno = va.getHighestImpactAnnotation();
		if (anno == null)
			return false; // no annotation for variant

		ImmutableSortedSet<VariantEffect> effects = anno.getEffects();
		String ann = anno.getAminoAcidHGVSDescription();
		if (effects.contains(VariantEffect.INTERGENIC_VARIANT)) {
			return true;
		} else if (effects.contains(VariantEffect.UPSTREAM_GENE_VARIANT)
				|| effects.contains(VariantEffect.DOWNSTREAM_GENE_VARIANT)) {
			// TODO(holtgrew): The annotation output format is different now, so this needs to be adjusted.
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
		} else if (effects.contains(VariantEffect.CODING_TRANSCRIPT_INTRON_VARIANT)
				|| effects.contains(VariantEffect.NON_CODING_TRANSCRIPT_INTRON_VARIANT)) {
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
	 * @param vc
	 *            VariantContext with the current annotation
	 * @param va
	 *            the annotations for <code>vc</code>
	 */
	private boolean doQualityControl(VariantContext vc, VariantAnnotations va) {
		if (this.doTargetFilter && isOffTarget(va)) {
			this.n_offTarget++;
			return false;
		}
		for (int i = 0; i < this.samples.size(); ++i) {
			SampleReport rep = this.samples.get(i);
			rep.addVariant(vc, va);
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
		JannovarData data = new JannovarDataSerializer(this.pathToSerializedJannovarFile).load();
		this.chromosomeMap = data.getChromosomes();
		this.refDict = data.getRefDict();
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