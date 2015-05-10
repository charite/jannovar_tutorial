# Minimal Maven Project Using Jannovar

## Get the Example Code

Perform the following steps to create an application that can use the jannovar-core library:

    $ git clone https://github.com/charite/jannovar_example_app.git
    $ cd minimal-jannovar-app.git
    $ mvn test
    ...

## Download a Jannovar Database

As a preparation, download the jannovar-cli application for downloading databases:

    $ wget https://github.com/charite/jannovar/releases/download/v0.14/jannovar-cli-0.14-bin.zip
    $ unzip jannovar-cli-0.14-bin.zip
    $ java -jar jannovar-cli-0.14/jannovar-cli-0.14.jar download hg19/refseq
    ...

## Build Example App

    $ mvn package
    ...
    $ java -jar java -jar target/example-project-0.1-SNAPSHOT.jar data/hg19_refseq.ser 
	Hello World!
	...

## JPed
JPed is a simple example application that is intended to help users get a feeling for how
realistic Jannovar-based applications can be written. The program takes as input a VCF file, a Ped file, 
and the transcript information file (Jannovar format), and filters all variants whose inheritance is
compatible with the mode of inheritance and the VCF file and the Ped file. After you have cloned this
repository (jannovar_tutorial), cd to the jped directory and enter

    $ mvn package
    $ java -jar target/jped-0.1-SNAPSHOT.jar 

This will show a usage message:

    Usage: java -jar JPed.jar -D ucsc_hg19.ser -V sample.vcf -P sample.ped  -I AR [-O fname]
    [INFO] where
    [INFO]
    [INFO] -D: the serialized transcript data file from Jannovar (e.g., ucsc_hg19.ser)
    [INFO] -V: the VCF file representing samples from a family
    [INFO] -P: the corresponding PED file
    [INFO] -I: the mode of inheritance to filter for (AD, AR, X, AR-HOM, AR-COMPHET)
    [INFO] -O: Name of output file (optional)

Consult the Jannovar manual for more information about the PED file format. Have a look
at the main function to understand how the code works.

* deserializeUCSCdata(): Input the transcript information file. This is the file created by the Jannovar application that contains information about each of the transcripts in the genome reference being used.
* annotateVCF(): Annotate each of the variants in the input VCF file
* processVariants(): extract one representative annotation for each variant (many of which have multiple annotations), and add it to a simple Gene object (the Gene will be the point of reference for deciding whether there is compatibility with the mode of inheritance; for instance, with autosomal recessive inheritance, if a gene has compound heterozygous variants in each affect individual, one of which is found in each parent, and no more than one of which is found in unaffected siblings, this would be compatible).
* parsePedFile(): obvious
* filterByInheritance(): Apply the rules encoded in the Jannovar library to decide if Genes/variants are compatible with the pedigree and the indicated mode of inheritance


## Have Fun!

From here on, it's your turn ;)