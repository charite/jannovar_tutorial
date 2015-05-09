package jqc;

import jannovar.common.Genotype;
import jannovar.common.VariantType;
import jannovar.exome.Variant;
import java.io.Writer;
import java.io.IOException;
import java.util.HashMap;


/**
 * Objects of this class encapsulate all of the data needed for
 * the Q/C report for a single sample in a VCF file.
 * @author Peter Robinson
 * @version 1.2 (December 29, 2013)
 */

public class SampleReport {
    /** ID or name of the sample in the VCF file. */
    private String sampleName=null;
    /** Index of this sample in the VCF file. */
    private int idx;

    /** Number of SNV variants in the VCF file. */
    private int n_SNV=0;
    /** Number of transversions in the VCF file. */
    private int n_transversions=0;
    /** Number of transitions in the VCF file. */
    private int n_transitions=0;
    /** Total number of called variants .*/
    private int n_totalVariants=0;
    /** Total number of synonymous single nucleotide variants. */
    private int n_synonymous=0;
    /** Total number of nonsynonymous single nucleotide variants. */
    private int n_nonsynonymous=0;
    /** Total number of homozygous variants. */
    private int n_homozygous=0;
    /** Total number of heterozygous variants. */
    private int n_heterozygous=0;
    /** Number of variants that are non-exoninc and non-splicing */
    private int n_offTarget=0;
    /** Number of variants called as "./.", i.e., missing */
    private int n_missing=0;
    /** Summary of the PHRED variant call score Q/C */
    private String phred=null;

    private String depth=null;
    /** Key: A variant type constant such as SPLICING.
     * (Constants are from jannovar.common.VariantType)
     * Value: A count for this variant type.
    */
    private HashMap<VariantType,Integer> variantTypeMap=null;
    /** A statistic object to calculate stats for the number of reads supporting a variant */
    private DescriptiveStatistic depthStat=null;

    /* Latex cannot deal with underscores, but these are common in sample names.
       Therefore, escape these characters.
    */
    private String escapeString(String s) {
	return s.replaceAll("_","\\\\_");
    }

    /**
     * Initializes the sample name and the variant type map.
     */
    public SampleReport(String name, int index) {
	this.sampleName = escapeString(name);
	this.idx = index;
	this.variantTypeMap = new HashMap<VariantType,Integer>();
	this.depthStat = new DescriptiveStatistic("Read depth");
    }

    /**
     * This method is intended to be called once
     * for every variant in the VCF file. It increments
     * the count of the respective variant type by one.
     */
    public void addVariant(Variant v) {
	Genotype gt = v.getGenotypeInIndividualN(this.idx); 

	if (v.isMissingInIndividualN(this.idx)) {
	    this.n_missing++;
	    return;
	}
	int rd = v.getVariantReadDepthIndividualN(this.idx);
	this.depthStat.addValue(rd);
	VariantType vt = v.getVariantTypeConstant();
	if (! this.variantTypeMap.containsKey(vt)) {
	    this.variantTypeMap.put(vt,1);
	} else {
	    Integer i = this.variantTypeMap.get(vt);
	    this.variantTypeMap.put(vt,i+1);
	}
	switch (gt) {
	case HOMOZYGOUS_REF:
	    /* Do nothing, not a variant in this individual */
	    break;
	case HOMOZYGOUS_ALT:
	    this.n_homozygous++;
	    break;
	case HETEROZYGOUS:
	    this.n_heterozygous++;
	    break;
	}
	if (v.isTransversion())
	    this.n_transversions++;
	else if (v.isTransition())
	    n_transitions++;
	if (v.isSynonymousVariant())
	    this.n_synonymous++;
	else if (v.is_missense_variant())
	    this.n_nonsynonymous++;
	if (v.isOffExomeTarget())
	    this.n_offTarget++;
	n_totalVariants++;
    }


    
    /** Adds the results of PHRED variant quality to this object. */
    public void setPhredString(String s) { this.phred = s;  }

    


    public void writeLatexSummary(Writer out) throws IOException {
	this.depthStat.evaluate();
	String offTarg = String.format("%d (%.1f\\%%)",n_offTarget, 100*(double)n_offTarget/n_totalVariants);
	out.write("\\section*{Sample: " + this.sampleName + "}\n\n");

	//out.write("A total of " + this.variantTypeMap.size() + " variants could be assigned to a variant class.\n");
	out.write("\\begin{center}\n");
	out.write("\\begin{tabular}{p{5cm}p{5cm}p{3cm}}\n");
	out.write("\\hline\n");
	out.write("Parameter & Value & Range\\\\ \n");
	out.write("\\hline\n");
	out.write("Variant call Phred score & " + this.phred +" & [$\\ >30$]\\\\ \n");
	out.write("Variant read depth & " + this.depthStat.toString() + " & -\\\\ \n");
	out.write("Total variants called & " + this.n_totalVariants + " & [$>25,000$]\\\\ \n");
	out.write("Transitions & " + this.n_transitions +"& - \\\\ \n");
	out.write("Transversions & " + this.n_transversions +"& - \\\\ \n");
	out.write(String.format("Ti/Tv  & %.3f & [2-3]  \\\\ \n", (double)this.n_transitions/this.n_transversions));
	if (this.n_missing > 0) {
	    out.write(String.format("Missing & %d (%.3f)& - \\\\ \n",this.n_missing,
				    100.0 * (double)n_missing / this.n_totalVariants));
	}
	out.write("Synonymous &  " + this.n_synonymous + " & -  \\\\ \n");
	out.write("Nonsynonymous &  " + this.n_nonsynonymous+ "& - \\\\ \n");
	out.write(String.format("ns/ss & %.2f & [0.8-1.0] \\\\ \n", (double)this.n_nonsynonymous/this.n_synonymous));
	out.write("Het calls & " + this.n_heterozygous + "& - \\\\ \n");
	out.write("Hom calls & " + this.n_homozygous + "& - \\\\ \n");
	out.write(String.format("Het/Hom & %.2f& [1.3-2.0] \\\\ \n", (double)this.n_heterozygous /this.n_homozygous)); 
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
    private void outputVariantRow(Writer out, VariantType vt) throws IOException {
	if (! this.variantTypeMap.containsKey(vt))
	    return;
	Integer c = this.variantTypeMap.get(vt);
	String var = VariantType.variantTypeAsString(vt);
	out.write(var + " & " + c + "\\\\ \n");
    }

	

	
	
    /**
     * Output a table of variant counts, classified accoding to 
     * variant type. */
    private void outputVariantDistribution(Writer out) throws IOException {

	out.write("\\begin{center}\n");
	out.write("\\begin{tabular}{p{6cm}p{4cm}}\n");
	out.write("\\hline\n");
	out.write("\\multicolumn{2}{l}{Variant Distribution}\\\\ \n");
	out.write("\\hline\n");
	out.write("Variant type & Count\\\\ \n");
	out.write("\\hline\n");
	outputVariantRow(out, VariantType.MISSENSE);
	outputVariantRow(out, VariantType.STOPGAIN);
	outputVariantRow(out, VariantType.FS_DELETION);
	outputVariantRow(out, VariantType.FS_INSERTION);
	outputVariantRow(out, VariantType.FS_SUBSTITUTION);
	outputVariantRow(out, VariantType.NON_FS_DELETION);
	outputVariantRow(out, VariantType.NON_FS_INSERTION);
	outputVariantRow(out, VariantType.NON_FS_SUBSTITUTION);
	outputVariantRow(out, VariantType.SPLICING);
	outputVariantRow(out, VariantType.STOPLOSS);
	outputVariantRow(out, VariantType.ncRNA_EXONIC);
	outputVariantRow(out, VariantType.ncRNA_SPLICING);
	outputVariantRow(out, VariantType.UTR3);
	outputVariantRow(out, VariantType.UTR5);
	outputVariantRow(out, VariantType.SYNONYMOUS);
	outputVariantRow(out, VariantType.INTRONIC);
	outputVariantRow(out, VariantType.ncRNA_INTRONIC);
	outputVariantRow(out, VariantType.UPSTREAM);
	outputVariantRow(out, VariantType.DOWNSTREAM);
	outputVariantRow(out, VariantType.INTERGENIC);
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