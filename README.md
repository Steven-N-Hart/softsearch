# SoftSearch

This is the new official repository to house SoftSearch, originally published in [PLoS One](http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0083356).  It is not under active development, but some people may find it useful.  To make deployments simple, we provided an installation script (`install.pl`), but the recommended way is using [Docker](https://www.docker.com/). 

## Option 1.  Run using DockerHub
```
docker run --rm -it stevenhart/softsearch perl softsearch/script/SoftSearch.pl [-cqlrmsd] -b <BAM> -f <Genome.fa>
	-q		Minimum mapping quality [20]
	-l		Minimum length of soft-clipped segment [5]
	-r		Minimum depth of soft-clipped reads at position [5]
	-m		Minimum number of discordant read pairs [5]
	-s		Number of sd away from mean to be considered discordant [6]
	-u		Number of unmated pairs [0]
	-d		Max distance between soft-clipped segments and discordant read pairs [Maximum normal insert]
	-o		Output file name [output.vcf]
	-t		Print temp files for debugging [no|yes]
	-c		use only this chrom or chr:pos1-pos2
	-p		use paired-end mode only. In other words, don't try to find soft-clipping events!
	-g		Enable paired-only seach. This will look for discordant read pairs even without soft clips.
	-a		set the minimum distance for a discordant read pair without soft-clipping info [10000]
	-L		Flag to print out even small deletions (low quality)
	-e		disable strict quality filtering of base qualities in soft-clipped reads [no]
	-blacklist	areas of the genome to skip calling.  Requires -genome_file
	-genome_file	tab seperated value of chromosome name and length.  Only used with -blacklist option

````

## Option 2. Create an environment to launch the container using Docker-Machine.

In this case, I'll just use `default`, but you can use whatever you want.
```
docker-machine create -d virtualbox default
eval $(docker-machine env default)
```
Next clone the repo from GitHub
```
git clone https://github.com/Steven-N-Hart/softsearch.git
cd softsearch
```
Now build the image
```
docker build -t test .

```
Finally, launch the container.  Note, that I am including a directory to mount.  This directory should include the paths to your `genome.fa` file and your `BAM` file
```
docker run --rm -it test perl softsearch/script/SoftSearch.pl [-cqlrmsd] -b <BAM> -f <Genome.fa>
	-q		Minimum mapping quality [20]
	-l		Minimum length of soft-clipped segment [5]
	-r		Minimum depth of soft-clipped reads at position [5]
	-m		Minimum number of discordant read pairs [5]
	-s		Number of sd away from mean to be considered discordant [6]
	-u		Number of unmated pairs [0]
	-d		Max distance between soft-clipped segments and discordant read pairs [Maximum normal insert]
	-o		Output file name [output.vcf]
	-t		Print temp files for debugging [no|yes]
	-c		use only this chrom or chr:pos1-pos2
	-p		use paired-end mode only. In other words, don't try to find soft-clipping events!
	-g		Enable paired-only seach. This will look for discordant read pairs even without soft clips.
	-a		set the minimum distance for a discordant read pair without soft-clipping info [10000]
	-L		Flag to print out even small deletions (low quality)
	-e		disable strict quality filtering of base qualities in soft-clipped reads [no]
	-blacklist	areas of the genome to skip calling.  Requires -genome_file
	-genome_file	tab seperated value of chromosome name and length.  Only used with -blacklist option

```

## Option 3.  Install on bare metal (not garunteed to work)
```
git clone https://github.com/Steven-N-Hart/softsearch.git
cd softsearch
perl install.pl -p $PWD

```