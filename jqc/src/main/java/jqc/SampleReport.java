package jqc;

import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;

import java.io.IOException;
import java.io.Writer;
import java.util.HashMap;

import com.google.common.collect.ImmutableSortedSet;

import de.charite.compbio.jannovar.annotation.Annotation;
import de.charite.compbio.jannovar.annotation.VariantAnnotations;
import de.charite.compbio.jannovar.annotation.VariantEffect;

/**
 * Objects of this class encapsulate all of the data needed for the Q/C report for a single sample in a VCF file.
 *
 * @author Peter Robinson
 * @version 1.2 (December 29, 2013)
 */

public class SampleReport {
	/** ID or name of the sample in the VCF file. */
	private String sampleName = null;
	/** Index of this sample in the VCF file. */
	private int idx;

	/** Number of SNV variants in the VCF file. */
	private int n_SNV = 0;
	/** Number of transversions in the VCF file. */
	private int n_transversions = 0;
	/** Number of transitions in the VCF file. */
	private int n_transitions = 0;
	/** Total number of called variants . */
	private int n_totalVariants = 0;
	/** Total number of synonymous single nucleotide variants. */
	private int n_synonymous = 0;
	/** Total number of nonsynonymous single nucleotide variants. */
	private int n_nonsynonymous = 0;
	/** Total number of homozygous variants. */
	private int n_homozygous = 0;
	/** Total number of heterozygous variants. */
	private int n_heterozygous = 0;
	/** Number of variants that are non-exoninc and non-splicing */
	private int n_offTarget = 0;
	/** Number of variants called as "./.", i.e., missing */
	private int n_missing = 0;
	/** Summary of the PHRED variant call score Q/C */
	private String phred = null;

	private String depth = null;
	/**
	 * Key: A variant type constant such as SPLICING. (Constants are from jannovar.common.VariantType) Value: A count
	 * for this variant type.
	 */
	private HashMap<VariantEffect, Integer> variantEffectMap = null;
	/** A statistic object to calculate stats for the number of reads supporting a variant */
	private DescriptiveStatistic depthStat = null;

	/*
	 * Latex cannot deal with underscores, but these are common in sample names. Therefore, escape these characters.
	 */
	private String escapeString(String s) {
		return s.replaceAll("_", "\\\\_");
	}

	/**
	 * Initializes the sample name and the variant type map.
	 */
	public SampleReport(String name, int index) {
		this.sampleName = escapeString(name);
		this.idx = index;
		this.variantEffectMap = new HashMap<VariantEffect, Integer>();
		this.depthStat = new DescriptiveStatistic("Read depth");
	}

	/**
	 * This method is intended to be called once for every variant in the VCF file. It increments the count of the
	 * respective variant type by one.
	 *
	 * @param vc
	 *            VariantContext with the current annotation
	 * @param va
	 *            the annotations for <code>vc</code>
	 */
	public void addVariant(VariantContext vc, VariantAnnotations va) {
		Genotype gt = vc.getGenotype(this.idx);
		Annotation anno = va.getHighestImpactAnnotation(); // only highest-impact one
		// tODO(holtgrew): Nick, is this behaviour desired here, I just jump out here in case of error
		if (anno == null)
			return; // no annotation

		if (gt.isNoCall()) {
			this.n_missing++;
			return;
		}
		// TODO(holtgrew): Nick, read depth will always be 0 here if not in VCF file
		int rd = (Integer) gt.getExtendedAttribute("AD", 0);
		this.depthStat.addValue(rd);
		ImmutableSortedSet<VariantEffect> ve = anno.getEffects();
		// TODO(holtgrew): Nick, the behaviour here changed to consider all sets in the effect set.
		for (VariantEffect eff : ve) {
			if (!this.variantEffectMap.containsKey(eff)) {
				this.variantEffectMap.put(eff, 1);
			} else {
				Integer i = this.variantEffectMap.get(eff);
				this.variantEffectMap.put(eff, i + 1);
			}
		}

		if (gt.isHomVar())
			this.n_homozygous++;
		else if (gt.isHet())
			this.n_heterozygous++;

		if ((gt.getAllele(0).getBaseString() == "A" && gt.getAllele(1).getBaseString() == "G")
				|| (gt.getAllele(0).getBaseString() == "G" && gt.getAllele(1).getBaseString() == "A")
				|| (gt.getAllele(0).getBaseString() == "C" && gt.getAllele(1).getBaseString() == "T")
				|| (gt.getAllele(0).getBaseString() == "T" && gt.getAllele(1).getBaseString() == "C"))
			this.n_transitions++;
		else if ((gt.getAllele(0).getBaseString() == "A" && gt.getAllele(1).getBaseString() == "C")
				|| (gt.getAllele(0).getBaseString() == "C" && gt.getAllele(1).getBaseString() == "A")
				|| (gt.getAllele(0).getBaseString() == "G" && gt.getAllele(1).getBaseString() == "T")
				|| (gt.getAllele(0).getBaseString() == "T" && gt.getAllele(1).getBaseString() == "G"))
			this.n_transversions++;

		if (ve.contains(VariantEffect.SYNONYMOUS_VARIANT))
			this.n_synonymous++;
		else if (ve.contains(VariantEffect.MISSENSE_VARIANT)) // TODO(holtgrew): Nick, this might be incorrect
			this.n_synonymous++;

		// TODO(holtgrew): Nick, I think you should check the implementatin of isOffExome to see whether it doe the right thing
		if (ve.first().isOffExome())
			this.n_offTarget++;

		n_totalVariants++;
	}

	/** Adds the results of PHRED variant quality to this object. */
	public void setPhredString(String s) {
		this.phred = s;
	}

	public void writeLatexSummary(Writer out) throws IOException {
		this.depthStat.evaluate();
		String offTarg = String.format("%d (%.1f\\%%)", n_offTarget, 100 * (double) n_offTarget / n_totalVariants);
		out.write("\\section*{Sample: " + this.sampleName + "}\n\n");

		//out.write("A total of " + this.variantTypeMap.size() + " variants could be assigned to a variant class.\n");
		out.write("\\begin{center}\n");
		out.write("\\begin{tabular}{p{5cm}p{5cm}p{3cm}}\n");
		out.write("\\hline\n");
		out.write("Parameter & Value & Range\\\\ \n");
		out.write("\\hline\n");
		out.write("Variant call Phred score & " + this.phred + " & [$\\ >30$]\\\\ \n");
		out.write("Variant read depth & " + this.depthStat.toString() + " & -\\\\ \n");
		out.write("Total variants called & " + this.n_totalVariants + " & [$>25,000$]\\\\ \n");
		out.write("Transitions & " + this.n_transitions + "& - \\\\ \n");
		out.write("Transversions & " + this.n_transversions + "& - \\\\ \n");
		out.write(String.format("Ti/Tv  & %.3f & [2-3]  \\\\ \n", (double) this.n_transitions / this.n_transversions));
		if (this.n_missing > 0) {
			out.write(String.format("Missing & %d (%.3f)& - \\\\ \n", this.n_missing, 100.0 * n_missing
					/ this.n_totalVariants));
		}
		out.write("Synonymous &  " + this.n_synonymous + " & -  \\\\ \n");
		out.write("Nonsynonymous &  " + this.n_nonsynonymous + "& - \\\\ \n");
		out.write(String.format("ns/ss & %.2f & [0.8-1.0] \\\\ \n", (double) this.n_nonsynonymous / this.n_synonymous));
		out.write("Het calls & " + this.n_heterozygous + "& - \\\\ \n");
		out.write("Hom calls & " + this.n_homozygous + "& - \\\\ \n");
		out.write(String.format("Het/Hom & %.2f& [1.3-2.0] \\\\ \n", (double) this.n_heterozygous / this.n_homozygous));
		out.write("Off target variants & " + offTarg + " & -  \\\\ \n");
		out.write("\\hline\n\n\n");
		out.write("\\end{tabular}\n");
		out.write("\\end{center}\n\n");
		out.write("\\vspace{1cm}\n\n");
		outputVariantDistribution(out);

	}

	/**
	 * Output one row of the variant count table (see {@link #outputVariantDistribution}).
	 */
	private void outputVariantRow(Writer out, VariantEffect ve) throws IOException {
		if (!this.variantEffectMap.containsKey(ve))
			return;
		Integer c = this.variantEffectMap.get(ve);
		String var = ve.getSequenceOntologyTerm().replace("_", "\\_");
		out.write(var + " & " + c + "\\\\ \n");
	}

	/**
	 * Output a table of variant counts, classified accoding to variant type.
	 */
	private void outputVariantDistribution(Writer out) throws IOException {

		out.write("\\begin{center}\n");
		out.write("\\begin{tabular}{p{6cm}p{4cm}}\n");
		out.write("\\hline\n");
		out.write("\\multicolumn{2}{l}{Variant Distribution}\\\\ \n");
		out.write("\\hline\n");
		out.write("Variant type & Count\\\\ \n");
		out.write("\\hline\n");

		// TODO(holtgrewe): Nick, you probably want to thin this out a bit... ;)
		outputVariantRow(out, VariantEffect.CODING_SEQUENCE_VARIANT);
		outputVariantRow(out, VariantEffect.CODING_TRANSCRIPT_INTRON_VARIANT);
		outputVariantRow(out, VariantEffect.CODING_TRANSCRIPT_VARIANT);
		outputVariantRow(out, VariantEffect.COMPLEX_SUBSTITUTION);
		outputVariantRow(out, VariantEffect.CONSERVED_INTERGENIC_VARIANT);
		outputVariantRow(out, VariantEffect.CONSERVED_INTRON_VARIANT);
		outputVariantRow(out, VariantEffect.DIRECT_TANDEM_DUPLICATION);
		outputVariantRow(out, VariantEffect.DISRUPTIVE_INFRAME_DELETION);
		outputVariantRow(out, VariantEffect.DISRUPTIVE_INFRAME_INSERTION);
		outputVariantRow(out, VariantEffect.DOWNSTREAM_GENE_VARIANT);
		outputVariantRow(out, VariantEffect.EXON_LOSS_VARIANT);
		outputVariantRow(out, VariantEffect.EXON_VARIANT);
		outputVariantRow(out, VariantEffect.FEATURE_TRUNCATION);
		outputVariantRow(out, VariantEffect.FIVE_PRIME_UTR_PREMATURE_START_CODON_GAIN_VARIANT);
		outputVariantRow(out, VariantEffect.FIVE_PRIME_UTR_TRUNCATION);
		outputVariantRow(out, VariantEffect.FIVE_PRIME_UTR_VARIANT);
		outputVariantRow(out, VariantEffect.FRAMESHIFT_ELONGATION);
		outputVariantRow(out, VariantEffect.FRAMESHIFT_TRUNCATION);
		outputVariantRow(out, VariantEffect.GENE_VARIANT);
		outputVariantRow(out, VariantEffect.INFRAME_DELETION);
		outputVariantRow(out, VariantEffect.INFRAME_INSERTION);
		outputVariantRow(out, VariantEffect.INITIATOR_CODON_VARIANT);
		outputVariantRow(out, VariantEffect.INTERGENIC_REGION);
		outputVariantRow(out, VariantEffect.INTERGENIC_VARIANT);
		outputVariantRow(out, VariantEffect.INTERNAL_FEATURE_ELONGATION);
		outputVariantRow(out, VariantEffect.INTRAGENIC_VARIANT);
		outputVariantRow(out, VariantEffect.INTRON_VARIANT);
		outputVariantRow(out, VariantEffect.MIRNA);
		outputVariantRow(out, VariantEffect.MISSENSE_VARIANT);
		outputVariantRow(out, VariantEffect.MNV);
		outputVariantRow(out, VariantEffect.NON_CODING_TRANSCRIPT_EXON_VARIANT);
		outputVariantRow(out, VariantEffect.NON_CODING_TRANSCRIPT_INTRON_VARIANT);
		outputVariantRow(out, VariantEffect.NON_CODING_TRANSCRIPT_VARIANT);
		outputVariantRow(out, VariantEffect.RARE_AMINO_ACID_VARIANT);
		outputVariantRow(out, VariantEffect.REGULATORY_REGION_VARIANT);
		outputVariantRow(out, VariantEffect.SEQUENCE_VARIANT);
		outputVariantRow(out, VariantEffect.SPLICE_ACCEPTOR_VARIANT);
		outputVariantRow(out, VariantEffect.SPLICE_DONOR_VARIANT);
		outputVariantRow(out, VariantEffect.SPLICE_REGION_VARIANT);
		outputVariantRow(out, VariantEffect.SPLICING_VARIANT);
		outputVariantRow(out, VariantEffect.START_LOST);
		outputVariantRow(out, VariantEffect.STOP_GAINED);
		outputVariantRow(out, VariantEffect.STOP_LOST);
		outputVariantRow(out, VariantEffect.STOP_RETAINED_VARIANT);
		outputVariantRow(out, VariantEffect.STRUCTURAL_VARIANT);
		outputVariantRow(out, VariantEffect.SYNONYMOUS_VARIANT);
		outputVariantRow(out, VariantEffect.TF_BINDING_SITE_VARIANT);
		outputVariantRow(out, VariantEffect.THREE_PRIME_UTR_TRUNCATION);
		outputVariantRow(out, VariantEffect.THREE_PRIME_UTR_VARIANT);
		outputVariantRow(out, VariantEffect.TRANSCRIPT_ABLATION);
		outputVariantRow(out, VariantEffect.TRANSCRIPT_VARIANT);
		outputVariantRow(out, VariantEffect.UPSTREAM_GENE_VARIANT);
		outputVariantRow(out, VariantEffect.FRAMESHIFT_VARIANT);

		out.write("\\hline\n");
		out.write("Total  & " + this.n_totalVariants + "\\\\ \n");
		out.write("\\hline\n");
		out.write("\\end{tabular}\n");
		out.write("\\end{center}\n\n");
		out.write("\\clearpage\n\n");
		out.write("\\newpage\n\n");
	}

}
/* eof */