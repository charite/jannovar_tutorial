package jped;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;

/** Command line parser from apache */
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.GnuParser;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.commons.cli.Parser;

import com.google.common.collect.ImmutableList;
import com.google.common.collect.ImmutableMap;

import de.charite.compbio.jannovar.JannovarException;
import de.charite.compbio.jannovar.annotation.Annotation;
import de.charite.compbio.jannovar.annotation.VariantAnnotations;
import de.charite.compbio.jannovar.annotation.VariantAnnotator;
import de.charite.compbio.jannovar.data.Chromosome;
import de.charite.compbio.jannovar.data.JannovarData;
import de.charite.compbio.jannovar.data.JannovarDataSerializer;
import de.charite.compbio.jannovar.data.ReferenceDictionary;
import de.charite.compbio.jannovar.htsjdk.InvalidCoordinatesException;
import de.charite.compbio.jannovar.htsjdk.VariantContextAnnotator;
import de.charite.compbio.jannovar.pedigree.Genotype;
import de.charite.compbio.jannovar.pedigree.GenotypeList;
import de.charite.compbio.jannovar.pedigree.GenotypeListBuilder;
import de.charite.compbio.jannovar.pedigree.PedFileContents;
import de.charite.compbio.jannovar.pedigree.PedFileReader;
import de.charite.compbio.jannovar.pedigree.Pedigree;
import de.charite.compbio.jannovar.pedigree.PedigreeExtractor;
import de.charite.compbio.jannovar.pedigree.Person;
import de.charite.compbio.jannovar.pedigree.compatibilitychecker.CompatibilityCheckerException;
import de.charite.compbio.jannovar.pedigree.compatibilitychecker.ad.CompatibilityCheckerAutosomalDominant;
import de.charite.compbio.jannovar.pedigree.compatibilitychecker.ar.CompatibilityCheckerAutosomalRecessive;
import de.charite.compbio.jannovar.pedigree.compatibilitychecker.ar.CompatibilityCheckerAutosomalRecessiveCompoundHet;
import de.charite.compbio.jannovar.pedigree.compatibilitychecker.xr.CompatibilityCheckerXRecessive;

/** Classes from the Jannovar library. */

/**
 * Example Jannovar program designed to demonstrate how to integrate simple pedigree filtering with exome variant
 * annotation.
 *
 * @author Peter Robinson
 * @version 1.2 (December 27, 2013)
 */
public class JPed {

	/** Map of Chromosomes (Jannovar object storing all transcripts from the Chromosome). */
	private ImmutableMap<Integer, Chromosome> chromosomeMap = null;
	/** Genome chromosome information */
	private ReferenceDictionary refDict = null;
	/** List of variants from input file to be analysed. */
	private ArrayList<VariantContext> variantList = null;
	/** List of all sample names of VCF file */
	private ArrayList<String> sampleNames = null;
	/** {@link VariantAnnotator} to use for generating variant annotations */
	private VariantContextAnnotator annotator = null;
	/**
	 * Pedigree of the persons whose samples were sequenced. Created on the basis of a PED file for multisample VCF
	 * files, or as a default single-sample Pedigree for single-sample VCF files.
	 */
	private Pedigree pedigree = null;
	/** Map of all Genes encountered during parsing of VCF file. The key is the gene symbol. */
	private HashMap<String, Gene> geneMap = null;
	/** Path to serialized Jannovar file with transcript information (e.g. ucsc.ser) */
	private String pathToSerializedJannovarFile = null;
	/** Path to a Variant Call Format (VCF) file that we will analyse. */
	private String pathToVCFfile = null;
	/** Path to a PED file (pedigree file) representing the samples on the VCF file. */
	private String pathToPedFile = null;
	/** Name of the outfile (default: jped.txt). */
	private String outfilename = "jped.txt";

	/**
	 * An enumeration to indicate the mode inheritance, with AR=autosomal recessive, AD=autosomal dominant, X=X
	 * chromosomal, and HOM=autosomal recessive (only homozygous), and COMPHET=autosomal recessive (only compound
	 * heterozygous).
	 */
	private enum inheritanceMode {
		AR, AD, X, HOM, COMPHET
	};

	/** Mode of inheritance */
	private inheritanceMode mode;

	/**
	 * An inner class desiged to hold all variants that correspond to a given gene.
	 */
	public class Gene {
		public String symbol;
		public ArrayList<Annotation> annots;
		public ArrayList<VariantContext> vars;

		public Gene(String s) {
			this.symbol = s;
			this.vars = new ArrayList<VariantContext>();
			this.annots = new ArrayList<Annotation>();
		}

		public void addVariant(VariantContext v, Annotation a) {
			this.vars.add(v);
			this.annots.add(a);
		}

		public String getAnnotations() {
			StringBuilder sb = new StringBuilder();
			sb.append(symbol + "\n");
			for (Annotation a : annots) {
				sb.append(a.getAminoAcidHGVSDescription() + "\n");
			}
			return sb.toString();
		}
	}

	public JPed(String argv[]) {
		parseCommandLineArguments(argv);
		this.geneMap = new HashMap<String, Gene>();
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
	 * This function is responsible for the input of the serialized transcript data file that is produced by the main
	 * Jannovar application and contains information about all transcripts The transcript models are put into
	 * {@link #chromosomeMap}, which will allow efficient searching.
	 */
	public void deserializeUCSCdata() throws JannovarException {
		JannovarData data = new JannovarDataSerializer(this.pathToSerializedJannovarFile).load();
		this.chromosomeMap = data.getChromosomes();
		this.refDict = data.getRefDict();
		this.annotator = new VariantContextAnnotator(data.getRefDict(), data.getChromosomes());
	}

	/**
	 * This function inputs a VCF file, and creates a list of all variants identified ({@link #variantList}) as well as
	 * of all samples names ({@link #sampleNames}).
	 */
	public void annotateVCF() throws JannovarException {
		VCFFileReader parser = new VCFFileReader(new File(this.pathToVCFfile), /* requireIndex= */false);
		this.variantList = new ArrayList<VariantContext>();
		for (VariantContext vc : parser)
			this.variantList.add(vc);
		parser.close();

		this.sampleNames = parser.getFileHeader().getSampleNamesInOrder();
	}

	/**
	 * For this example, we group the variants according to the gene they correspond to. In a real program, you would
	 * implement the logic of your analysis in a function(s) like this, and filter variants according to your criteria
	 * of interest. For instance, you might remove variants that are common in the population by comparing them with
	 * data from dbSNP, or you might prioritize the variants by checking their function with GO data or comparing the
	 * phenotype of your patient with HPO data.
	 *
	 * @throws InvalidCoordinatesException
	 *             on problems with coordinates
	 */
	public void processVariants() throws JannovarException, InvalidCoordinatesException {
		int c = 0;
		for (VariantContext vc : this.variantList) {
			for (VariantAnnotations annos : annotator.buildAnnotations(vc)) {
				Annotation anno = annos.getHighestImpactAnnotation();
				if (anno == null)
					continue; // skip empty annotation
				if (anno.getEffects().first().isOffExome())
					continue; // skip intronic and intergenic
				String sym = anno.getGeneSymbol();
				Gene g = null;
				if (this.geneMap.containsKey(sym)) {
					g = this.geneMap.get(sym);
				} else {
					c++;
					if (c % 500 == 0) {
						System.out.println("Gene " + c + ": " + sym);
					}
					g = new Gene(sym);
					this.geneMap.put(sym, g);
				}
				g.addVariant(vc, anno);
			}
		}
	}

	/**
	 * Inputs a Ped file for multi-sample VCF files. If necessary, also adjust the order of samples in the Pedigree to
	 * match the order in the VCF file.
	 *
	 * @throws IOException
	 *             in the case of I/O errors
	 */
	public void parsePedFile() throws JannovarException, IOException {
		PedFileContents contents = new PedFileReader(new File(this.pathToPedFile)).read();
		String pedName = contents.getIndividuals().get(0).getPedigree();
		ImmutableList<Person> members = new PedigreeExtractor(pedName, contents).run();
		this.pedigree = new Pedigree(pedName, members);//.subsetOfMembers(this.sampleNames);
	}

	public void filterByInheritance() throws CompatibilityCheckerException {
		try {
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

			out.write(this.pedigree.toString() + "\n");

			for (Gene g : this.geneMap.values()) {
				GenotypeList lst = constructGenotypeList(g.vars, g.symbol);
				if (this.mode == inheritanceMode.AR && new CompatibilityCheckerAutosomalRecessive(pedigree, lst).run()) {
					out.write(g.getAnnotations());
				} else if (this.mode == inheritanceMode.COMPHET
						&& new CompatibilityCheckerAutosomalRecessiveCompoundHet(pedigree, lst).run()) {
					out.write(g.getAnnotations());
				} else if (this.mode == inheritanceMode.HOM
						&& new CompatibilityCheckerAutosomalRecessive(pedigree, lst).run()) {
					out.write(g.getAnnotations());
				} else if (this.mode == inheritanceMode.AD
						&& new CompatibilityCheckerAutosomalDominant(pedigree, lst).run()) {
					out.write(g.getAnnotations());
				} else if (this.mode == inheritanceMode.X && new CompatibilityCheckerXRecessive(pedigree, lst).run()) {
					out.write(g.getAnnotations());
				}
			}
			out.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	private GenotypeList constructGenotypeList(ArrayList<VariantContext> vars, String geneSymbol) {
		GenotypeListBuilder builder = new GenotypeListBuilder(geneSymbol, vars.get(0).getSampleNamesOrderedByName());
		String chrName = vars.get(0).getChr();
		builder.setIsXChromosomal(refDict.getContigNameToID().containsKey("X")
				&& refDict.getContigNameToID().get(chrName) == refDict.getContigNameToID().get("X"));

		for (VariantContext vc : vars) {
			ImmutableList.Builder<Genotype> gts = new ImmutableList.Builder<Genotype>();
			for (String sample : vars.get(0).getSampleNamesOrderedByName()) {
				htsjdk.variant.variantcontext.Genotype genotype = vc.getGenotype(sample);
				if (genotype.getAllele(0).equals(vc.getReference()) && genotype.getAllele(1).equals(vc.getReference()))
					gts.add(Genotype.HOMOZYGOUS_REF);
				else if (!genotype.getAllele(0).equals(vc.getReference())
						&& !genotype.getAllele(1).equals(vc.getReference()))
					gts.add(Genotype.HETEROZYGOUS);
				else
					gts.add(Genotype.HOMOZYGOUS_REF);
			}
			builder.addGenotypes(gts.build());
		}

		return builder.build();
	}

	/**
	 * Parse the command line using Apache's CLI. A copy of the library is included in the Jannovar tutorial archive in
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
			options.addOption(new Option("P", "ped", true, "Path to ped file"));
			options.addOption(new Option("I", "inheritance", true, "Inheritance pattern (AR, AD, X)"));
			options.addOption(new Option("O", "outfile", true, "name of outfile (default: jlink.txt)"));
			Parser parser = new GnuParser();
			CommandLine cmd = parser.parse(options, args);
			if (cmd.hasOption("H")) {
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
	 * Parse a command line argument representing the inheritance mode (should be one of AR, AD, or X; otherwise the
	 * program will terminate with an error).
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
			System.out
					.println("Please use one of \"AR\", \"AR-HOM\", \"AR-COMPHET\", \"AD\", or \"X\", e.g., \"-I AD\"");
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
