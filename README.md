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

## Have Fun!

From here on, it's your turn ;)