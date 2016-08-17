# README #

This README would normally document whatever steps are necessary to get your application up and running.

### What is this repository for? ###

This repository exits to make RMAExtractor available to interested and/ or contributing parties

### How do I get set up? ###

You have to choose a directory to store the ncbi.map and ncbi.tre file located in the 
resources folder
You than have to change the path variables in of NCBI_MapReader and NCBI_TreeReader 
currently located in the NCBI_MapReader folder to point to this location. 
Lastly you have to export the project as executable jar and give read and write permissions
to the jar.
the program should be executed via the folowing steps
Also use the mutithreaded_alpha branch

java  -jar nameOfJar.jar Path/To/Directory/containing/RMA6.files 
Path/to/output/directory/ 
Path/to/file/taxonlist (has to end with taxon.txt) 
topPercent(value as double 0.1 not .1 or 1) 
Number_of_CoresAsInteger (eg 4 by default 1 at maximum all available to java runtime)

it is highly advised to specify an output dir otherwise the output files will be written into the input directory!!!
The summary for all files will be called overallSummary.txt
supplemetry information will start with the filename+readDist.text or supplement.txt.gz
### Contribution guidelines ###

* Writing tests
* Code review
* Other guidelines

### Who do I talk to? ###

* Repo owner or admin
* Other community or team contact