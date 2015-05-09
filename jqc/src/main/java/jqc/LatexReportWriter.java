package jqc;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.Writer;

/**
 * A simple class to coordinate writing results of QC to a latex file.
 *
 * @author Peter Robinson
 * @version 0.2 (December 31, 2013).
 */
public class LatexReportWriter {

	private String vcfFileName = null;

	public LatexReportWriter(String vcfName) {
		this.vcfFileName = vcfName;
	}

	/**
	 * Try to run pdflatex from the command line to generate a PDF file automatically.
	 */
	public void createPDFFile(String path) {
		Runtime runtime = Runtime.getRuntime();
		String cmd = String.format("pdflatex %s", path);
		System.out.println("[INFO] Attempting to create PDF file with command\n[INFO]  *--:" + cmd + "--*");
		System.out.println("[INFO] If this commandhangs, automatic generation of PDF file did not work. ");
		System.out.println("[INFO] Stop operation by ctl-z and compile by hand with pdflatex");
		try {
			Process p = Runtime.getRuntime().exec(cmd);
			BufferedReader bre = new BufferedReader(new InputStreamReader(p.getErrorStream()));
			String line = null;
			while ((line = bre.readLine()) != null) {
				System.err.println(line);
			}
			bre.close();
			p.waitFor();
			//int exitval = p.exitValue();
			System.out.println("[INFO]  Done generating PDF output.");
		} catch (Exception e) {
			System.err.println("Error encountered while trying to generate PDF file");
			System.err.println("JQCVCF generated a latex file that can be compiled to a PDF file with pdflatex");
			e.printStackTrace();
			System.exit(1);
		}

	}

	/**
	 * Write the header and some introductory information to a latex file.
	 */
	public void writeHeader(Writer out) throws IOException {
		String fname = this.vcfFileName.replaceAll("_", "\\\\_");
		out.write("\\documentclass[10pt]{article}\n" + "\\usepackage[utf8]{inputenc}\n"
				+ "\\usepackage[english]{babel}\n" + "\\usepackage{fontenc}\n" + "\\usepackage[margin=2cm]{geometry}\n"
				+ "\\title{Exome Q/C Report" + fname + "}\n" + "\\begin{document}\n" + "\\section*{Exome Q/C Report: "
				+ fname + "}\n\n");

	}

	/**
	 * Write a description of the overall quality control results. (latex subsection)
	 */
	public void description(Writer out, int n_samples, int n_totalVariants) throws IOException {
		out.write("\n\n\n");
		out.write("\\subsection{Total VCF data}\n");
		if (n_samples > 1) {
			out.write("There were " + n_samples + " samples and a total of " + n_totalVariants
					+ " identified variants.\n");
		} else {
			out.write("There were a total of " + n_totalVariants + " identified variants.\n");
		}
	}

	public void outputOffTargetSummary(Writer out, int distanceThreshold, int n_offTarget, int n_total)
			throws IOException {
		String off = String.format("%.1f", 100 * (double) n_offTarget / n_total);
		out.write("\\begin{itemize} \\item Variants were filtered to be with " + distanceThreshold
				+ " nucleotides of an exon. " + n_offTarget + " off-target variants (" + off + "\\%) were "
				+ " removed before q/c analysis. \\end{itemize}\n");
	}

	public void writeExplanation(Writer out) throws IOException {
		out.write("\\paragraph{Explanations}\n");
		out.write("\\begin{scriptsize}\n");
		out.write("\\paragraph{Phred Score}\n");
		out.write("By convention, a good Phred score is taken to be anything above 30 for base calls. For variant calls, the actual scores may be much higher and the standard range differs from program to program and platform to platform.\n");
		out.write("\\paragraph{Total variants called}\n");
		out.write("This depends highly on sequencing depth, but is often above 30,000 total variants\n");
		out.write("\\paragraph{Ti/Tv}\n");
		out.write("Ti/Tv describes the transition to transversion ratio, and should be roughly in the range of 	2--3. ");
		out.write("Note that there are twice as many possible transversions as transitions, and therefore, if there were no biological bias one might expect Ti/Tv=0.5. Mark DePristo states that Ti/Tv should be 0.5 for false positives i.e. under randomness / absent biological forces, implying that it is not adjusted but rather simply the ratio of the number of transitions to the number of transversions. In practice, the Ti/Tv ratio should be about 2.1 for whole genome sequencing and 2.8 for whole exome, and that if it is lower, it may suggest the data includes false positives caused by random sequencing errors. The greater transition content of the exome results from the fact that the exome is under stronger selective pressure against missense mutations, whereas many transitions are tolerated as silent mutations because transitions in the third base of a codon rarely change the amino acid.\n");
		out.write("\\paragraph{ns/ss}\n");
		out.write("The ratio of non-synonymous substitutions to synonymous substitutions. Since synonymous substitutions are better tolerated by evolution,  the ratio may be expected to be below 1  but random errors would skew it upward. There are no accepted normal values, but we have generally seen values between 0.8 and 1.0.\n");
		out.write("\\paragraph{het/hom}:\n");
		out.write("The ratio of heterozygous to homozygous variant genotypes across all sample-site combinations. There are no accepted normal values, but populations with recent admixture will skew towards heterozygosity, whereas populations with inbreeding will skew towards homozygosity. We have generallyseen values between 1.3 and 2.0.\n");
		out.write("\\paragraph{missing}:\n   This refers to the number of variants called as \"./.\". These are variants that were called in a sample of a multisample VCF file but not in the current sample.\n");
		out.write("\\paragraph{References}\n");
		out.write("\\begin{itemize}\n");
		out.write("\\item Bainbridge MN et al. (2011) Targeted enrichment beyond the consensus coding DNA sequence exome reveals exons with higher variant densities. {\\it Genome Biol}. {\\bf 12}:R68\n");
		out.write("\\item Tennessen JA et al (2012) Evolution and functional impact of rare coding variation from deep sequencing of human exomes. {\\it Science}. {\\bf 337}:64-9.\n");
		out.write("\\item J\\\"{a}ger M et al (2014) Jannovar: A Java library for Exome Annotation.\n");
		out.write("\\end{itemize}\n");
		out.write("\\end{scriptsize}\n");
		out.write("\\paragraph{About}\n");
		out.write("This Quality Control report was generated " + " using JQCVCF from the Jannovar project. ");
		out.write("See the following websites for more information: \\newline \n"
				+ "https://github.com/charite/jannovar and \\newline \n"
				+ "http://compbio.charite.de/contao/index.php/jannovar.html.\n");
	}

	/**
	 * Write the footer (final line) of the latex document.
	 */
	public void writeFooter(Writer out) throws IOException {
		writeExplanation(out);
		out.write("\\end{document}\n");
	}

}