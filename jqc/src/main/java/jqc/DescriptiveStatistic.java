package jqc;


import java.util.ArrayList;
import java.util.Collections;


/**
 * This is a convenience class that we will use to calculate the mean, median, and standard deviation
 * of some of the qualities we are using for Q/C of VCF files.
 * @author Peter Robinson
 * @version 1.0 (November 1, 2013)
 */

public class DescriptiveStatistic {

    private ArrayList<Double> values=null;

    private double mean;

    private double median;

    private double variance;

    private double standardDeviation;

    private String label=null;
    
    /**
     * @param lab A label, e.g., "PHRED Score".
     */
    public DescriptiveStatistic(String lab) {
	this.label = lab;
	this.values = new ArrayList<Double>();
    }

    public void addValue(double v) {
	this.values.add(v);
    }

    /**
     * This function should be called after all values have been entered.
     * It calculates the mean, median, and standard deviation.
     */
    public void evaluate() {
	Collections.sort(this.values);
	calculateMeanAndSD();
	calculateMedian();
    }


    private void calculateMeanAndSD() {
	if (this.values.size()==0) {
	    System.err.println("Warning attempt to calculate mean for zero-size list");
	    System.exit(1);
	}
	double sum = 0;
	double squared_sum = 0;
	int N = this.values.size();
	
	for (int i = 0; i < N; i++) {
	    sum += this.values.get(i);
	    squared_sum += this.values.get(i) * this.values.get(i);
	}
	this.mean = sum / this.values.size();
	this.variance = squared_sum / N - mean * mean;
	this.standardDeviation = Math.sqrt(this.variance);
    }

    /** This function assumes that 
     * {@link #values} has been sorted.
     */
    private void calculateMedian() {
	Collections.sort(this.values);
	int len = this.values.size();
	int middle = len/2;
	if (len%2 == 1) {
	    this.median = this.values.get(middle);
	} else {
	    this.median= (this.values.get(middle-1) + this.values.get(middle)) / 2.0;
	}
    }


   

    public String toString() {
	String s = String.format("$%.1f\\pm %.1f$; median: %.1f",
				 this.mean,this.standardDeviation,this.median);
	return s;
    }


}