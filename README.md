# README #
README

This README is actually quite outdated. Please check the HOPS readme for current information

How to run MaltExtract:

```
java –Xmx20G –jar MaltExtract.jar 
-i(nput) /myExampleRMAs/
-o(utput)/myOutDirectory/
-f(ilter) default
-t(axa) “Yersinia pestis”
-p(threads) 20
```

MaltExtract is a Java program and dependent on Java version 1.8 or higher.
So you have to communicate to your Operating system (Unix or OSX) that you want to run a java program you do that but typing java into the command line. Next you have the option set java up with special parameters we will for example increase the maximum amount of memory java can allocate with –xmX20G to 20 Gigabytes. Next we tell Java which type of program we want to run. In our case a jar with –jar. Next we provide the name/path of the MaltExtract.jar.

Now we can provide MaltExtract with its parameters.
For every run it is necessary to provide an input file or directory. We can do this with –i (or –input, if you like to write more) followed by the name of the RMA6 file or directory containing the files. In our case we write –i /myExampleRMAs/.
Additionally you always have to provide a directory for the output. This is done with –o (or output) and the name of the directory.  The directory path has to follow Operating system standards. You are not allowed to use C:\myOutput\ in a Unix environment. We write -o /myOutDirectory/. MaltExtract will create the directory for you if it does not already exit.
Next we have to specify which filter we want to use.  This is done with –f (or filter) and the name of the chosen filter. We choose default (which does not do any filtering but does produce additional proofs, like edit distance distribution or read length distribution). We therefore write –f default. Since we use a filter that actually has to open the file and go through it we also have to specify a list of target species.  We choose Yersinia pestis with –t “Yersinia pestis”.  
MaltExtract is multithreaded so we can speed up the analysis by assigning more cores to the program. With –p 20 we assign 20 cores to MaltExtract to allow for one file per core.
Now it is time to grab a nice cup of coffee (or can I recommend a good cup of tea) and let MaltExtract work its magic. When it’s done we will talk about how to plot part of the output and find which samples are actually significant. 


# Manual :

For using MaltExtract in a multithreaded fashion increases the amount of memory it needs to work efficient. So it is a good idea to increase the heap size of java with –xmX(NumberOfCores)G, (e,g –xmX50G). However the memory requirments vary based on sequencing depth and the amount of Reads that are assigned to species in the tagrget list as well as the amount of files that are procesed simultaniously.If the RMA6 files are very big it might be advisable to increase heap size even more. MaltExtract will report the amount of consumed mamor ar the end.
MaltExtract’s major parameters generally have a short and a long name. More advanced parameters generally only have a long version.
Running in MaltExtract in with out additional parameters except input and output will produce a summery tab delaminated summary file that lists that assigned read count for each node in each input file. This is set as default behavior, as it is the quickest filter and it is useful information to have for any further analysis. At the moment the idea behind MaltExtract to use MALT in a very lenient fahion to have an indication on Backgorund noise and do stricter filtering with RMAE

-i(nput) required parameter to specify the input files, directories or symbolic links to files(in which case the name of the link will be used instead of the source file name). At least one file or directory has to be provided in order to start MaltExtract

-o(utput) required parameter. An output directory has to be provided in order to use MaltExtract. 

-t(axa)  provide either one  or many Megan species name(s) in “”. Like “Yersinia pestis” “Bordetella” “Mycobacterium tuberculosis” as a white space separated line. Or provide the path to a file containing a new line separated list of species.
This is a required parameter for every filter except –scan.

-f(ilter) filter is used to trigger preset behaviors that affect which reads are considered for all outputs.  “
default” does not do any actual filtering but produces additional outputs in addition to the runSummary.txt. The output will be saved in a subfolder in the output directory named after the filter.

“ancient” considers only reads which have at least one alignment in the top percent of alignments that either shows a C>T substitution in the first 5 positions from the 5’ end or a G>A in the first 5 positions of the 3’ end. These substitutions are very common in ancient DNA. 

“def_anc” (or "full") simultaneously runs the default and the ancient filter

“scan” quickly produces an unfiltered summary file into the output directory does not produce additional outputs. Is the only option that works without setting the parameter taxa to allow to retrieve all assigned nodes. 

"crawl" retrieve for a species that damage pattern, edit distance, percent identity and read length distribution of all strains/alignments of the -t Species flag - so don’t need to re-map with BWA for mapDamage. Now this process is generally time intensive and should only used scarcely on a smaller subset of interesting files

-p(threads) MaltExtract is multithreaded in the sense that different files can be analyzed on different cores. Assigning more cores to MaltExtract can speed up the analysis of more than one file. However assigning even more cores than there are files to be processed will not further increase the speed. 

-resources can be used to specify a directories containing another ncbi.tre and ncbi.map file These files contain the NCBI’s phylogeny and the assigned Megan IDs for each species. They will be frequently updated in the Daniel Houstan’s Megan Megan github page under resources. However in the current version of MaltExtract this parameter is no longer required as MAltExtract will simply download the necessary files from my github

-top use top 0.XX percent of alignments, default uses value used in MALT, e.g. 10% best alignments maintained. However it is possible to be either stricter or more lenient than that.

-maxReadLength use only reads shorter than X - note, all reads out of ClipAndMerge/EAGER already >30bp. This option is disabled by default. This can be used with ancient and reads to output reads that could satisfy the conditions to be counted as ancient (short length C>T and G>A substitutions depending on UDG treatment and sequencing method)

-minPI only consider Alignments with a min percent identity XX.XX. This option is disabled by default as well and at the moment the minimum percent identity actually present in RMA6 files is 85.00. However this allows users to Malt data very unspecific and be stricter when using MaltExtract.

-minComp only use reads with Complexity higher as 0.XX, Wootton & Federhen complexity, string of monobase == 0, 50% GC/AT == 0.8 (correlates to Shannon entropy). Can be used to ignore low complexity reads that map very unspecific.  A useful value is around 0.6- 0.75. However this filter is normally disabled. Well to be precise the function calculates the mean complexity for 16bp sliding windows. 

-reads output alignments in BLAST format to a file in in the out directory’s subfolder reads. Only encouraged for nodes with a small number of reads. Disabled by default.

--destackingOff  turn off stacked read removal which is usually turned on for all filters and removes alignments on a reference that overlap at least one position.
Filter will automatically turn itself of if average coverage on a reference is higher than 10

--dupRemOff turn off duplicate removal. Normally if reads share the exact same start and end position as another read only the first read is kept. 

--downSampOff turn off downsampling to 10000 per node. Turining off downsampling will slightly increase runtime but will also increase sensitivity

--useTopAlignment Only use the highest Scoring Alignment per Read will usually improve Sensitivity and somewhat decrease runtime. 

--meganSummary request Megan summury files that are readable with MEGAN and  contain the entries of the rma6 file and the unfiltered number of reads assigned to a node

-v(erbose) if chosen MaltExtract will print more detailed output to the screen and its logs

-h(elp) print all parameters, accepted values and usage to the screen than orderly shut down the program. 


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

```
java  -jar nameOfJar.jar Path/To/Directory/containing/RMA6.files 
Path/to/output/directory/ 
Path/to/file/taxonlist (has to end with taxon.txt) 
topPercent(value as double 0.1 not .1 or 1) 
Number_of_CoresAsInteger (eg 4 by default 1 at maximum all available to java runtime)
```

it is highly advised to specify an output dir otherwise the output files will be written into the input directory!!!
The summary for all files will be called RunSummary.txt
supplemetry information will start with the filename+AlignmentDist.text or supplement.txt.gz


### Explanation of Output Files sorted by Folder ###

For read distribution and coverage ane reference sequence is chosen to be representive. This will always be the reference with the most alignments.
This is necessary as Malt allows reads to have alignments to mutliple references. This however is only a problem for nodes that are not in the leaves of the taxonomic tree.
there is a addtional node entries file in the alignment folders to provide some indication on how other files behaved on the same node

If you use --useTopALignment flag in MaltExtract than all output except read distribution will reflect only all topalignments of these noodes,
that option is currently the default in MAltExtract1.4


# coverage #

_coverageHist.txt
Taxon	Reference	0	1	2	3	4	5	6	7	8	9	10	higher
Actinomyces_oris	Actinomyces_oris	1465520	311419	310106	251124	194992	138740	99276	69300	54112	41976	27670	78682

this file is available for each input file and starts by telling the taxonomic node we are analyzing, than the reference that has the most assinged alignments. Lastly we have the number of positions
covered 0,1 ,2 .... and so on times while anything covered higher than 10 times is summarized in the higher bin. As HOPS is desinged as a screening tool and not for higher covered data. Spikes in
coverage are usually a sign of something going wrong. FOr example reads mappng to conserved regions due to reference bias

_postionsCovered.txt higher than X

Taxon	Reference	AverageCoverage	Coverge_StandardDeviation	percCoveredHigher1	percCoveredHigher2	percCoveredHigher3	percCoveredHigher4	percCoveredHigher5
Actinomyces_oris	Actinomyces_oris	0.518	2.027	0.221	0.119	0.068	0.04	0.024


# damage mismatch #

_damageMismatch.txt
Node	C>T_1	C>T_2	C>T_3	C>T_4	C>T_5	C>T_6	C>T_7	C>T_8	C>T_9	C>T_10	G>A_11	G>A_12	G>A_13	G>A_14	G>A_15	G>A_16	G>A_17	G>A_18	G>A_19	G>A_20	D>V(11Substitution)_1	D>V(11Substitution)_2	D>V(11Substitution)_3	D>V(11Substitution)_4	D>V(11Substitution)_5	D>V(11Substitution)_6	D>V(11Substitution)_7	D>V(11Substitution)_8	D>V(11Substitution)_9	D>V(11Substitution)_10	H>B(11Substitution)_11	H>B(11Substitution)_12	H>B(11Substitution)_13	H>B(11Substitution)_14	H>B(11Substitution)_15	H>B(11Substitution)_16	H>B(11Substitution)_17	H>B(11Substitution)_18	H>B(11Substitution)_19	H>B(11Substitution)_20	considered_Matches
Actinomyces_oris	0.38113147410358567	0.22811386778474826	0.1265068493150685	0.10326872039506836	0.08713955623081063	0.07130964610861171	0.05629024440905629	0.05195253641386238	0.05465227490094275	0.04395488932474082	0.0512072707542051	0.05394391490537367	0.05460408904474154	0.05873861488228218	0.0754444640028164	0.08394833948339483	0.10228668941979523	0.12681842672413793	0.22373607274758697	0.3721667898497167	0.017623124047259124	0.013325395370910828	0.009091226891895756	0.007938668646919925	0.008321089167541502	0.007666419135009274	0.006591192463455534	0.007089080675899139	0.006965138291099348	0.005887792946301167	0.005970421202834361	0.006467250078655744	0.007429127731631897	0.006685473422833154	0.007647351075809306	0.008448209562207954	0.007967270735719878	0.009546741639450544	0.013146367481755575	0.017408078712948378	85817

The first value is taxonomic unit
the next 10 values are C>T substitutions in 5' direction for all alignments for all reads
the next 10 values are G>A substitutions in 3' direction for all alignments for all reads
the next 10 values are Non C>T substitutions in 5' direction for all alignments for all reads which allows us to estimate noise
the next 10 values are Non G>A substitutions in 3' direction for all alignments for all reads which allows us to estimate noise
the last node is the number of processed alignments

# editDistance #

_editDistance.txt
Taxon	0	1	2	3	4	5	6	7	8	9	10	higher
Actinomyces_oris	3351	7926	11330	12708	11302	10548	8675	6875	4642	3350	2076	3974

First column is the target node, returns the editDistance Distribution for all topalingments of reads that are kept after filtering


# filterInformation #

_filterTable.txt
this file gives some summary statistics for all filtered reads and alignments.
First column is the taxonomic unit
Node	NumberOfUnfilteredReads	NumberOfFilteredReads	NumberOfUnfilteredAlignments	numberOfAlignments	turnedOn?
Actinomyces_oris	520523	520523	568637	520523	On

the next two columns are the number of reads prior to filtering and after filtering the difference between those numbers is the number of removed reads
the next two columnsare the number of alignments prior and after filtering, the difference between those aolumns is the number if removed alignments
the last column gives information on whether destacking was turned on or not (moght be bugged at the moment)


# percentIdentity #

_percentIdentity.txt

Taxon	80	85	90	95	100
Actinomyces_oris	0	13569	36455	31446	5287
First Column is the taxonomic Unit and  the next  5
columns gives the number of top alginments of filtered reads that fall into the set percent identitiy bins


# readDist #

_additionalNodeEntries.txt
TargetNode	01	02	03	04	05	06	07	08	09	10
Actinomyces_oris	Actinomyces_oris;_CP014232.1;_TOPREFPERCREADS100	Actinomyces_naeslundii;_AP017894.1;_TOPREFPERCREADS005	Actinomyces_oris_K20;_AB573870.1;_TOPREFPERCREADS002	Actinomyces_naeslundii;_EU621354.1;_TOPREFPERCREADS001	Actinomyces_naeslundii;_EU621259.1;_TOPREFPERCREADS001	Actinomyces_viscosus;_EU621357.1;_TOPREFPERCREADS001	Actinomyces_viscosus;_EU620893.1;_TOPREFPERCREADS001	Actinomyces_sp._Chiba101;_AP017896.1;_TOPREFPERCREADS001	Actinomyces_viscosus;_EU621241.1;_TOPREFPERCREADS000	Actinomyces_johnsonii;_EU621007.1;_TOPREFPERCREADS000

This File exits to provide a way for the user to infer how reads distirbute across multiple reference. Reads assigend higher than strain level will have alignments to multiple reference
the column are the 10 refernes with the most assinged alignments if the taxonomic unit has less than 10 references sequences the mossing entries will be replaced with a NA
Each column gives the name of the reference, and the perceantage of alignments assingend to it (when compared to the highest scoring reference)

_alignmentDist.txt
Taxon	Reference	uniquePerReference	nonStacked	nonDuplicatesonReference	TotalAlignmentsOnReference	ReferenceLength
Actinomyces_oris	Actinomyces_oris	0.197	3298	26476	86455	3042917
This file provides information on alignment distribution for the reference sequence that has the most alignments for each taxonomic unit.
The first two columns is the the name of the taxon and the name of the referece squence
the next is score that gives the fraction of bp that are covered uniquely over all the basepairs of the reads that could theoretically be uniquely mapping.
So this score gives some change to estimate overlap. the next three columns is the number of non stacking reads, of non duplicate reads and the total number of alignments assigned to the reference.
the last one is the lenght of the reference

Taxon	25bp	30bp	35bp	40bp	45bp	50bp	55bp	60bp	65bp	70bp	75bp	80bp	85bp	90bp	95bp	100bp	105bp	110bp	115bp	120bp	125bp	130bp	135bp	140bp	145bp	150bp	155bp	160bp	165bp	170bp	175bp	180bp	185bp	190bp	195bp	200bp
Actinomyces_oris	0	1917	7111	10570	11733	10835	8965	7580	6250	4657	4811	3022	2273	1835	1117	962	716	655	548	407	277	245	177	94	0	0	0	0	0	0	0	0	0	0	0	0
the last entriy provides for all filtered assinged reads how read length was distributed between 20 and 200 bp 

# RunSummary.txt #

This file returens for all input files how many reads are assinged to the target taxonomic units after filtering

# TotalCount.txt #

contains the total number of assinged reads for all files

### Who do I talk to? ###
