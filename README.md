# LycSpeciationGenomics
Notes and code for genomic and phenotypic analyses of the speciation continuum in Lycaeides butterflies. This is a fork of the contemporary evolution data set and analyses that excludes the time dimension.

 See [Pool_DNA_plates](https://drive.google.com/drive/folders/1U4AsshyMvlySNtODuSLWo0dYDH_rgho4) for details on samples and sample concentrations.

# Data
The pool-seq data are currently in /uufs/chpc.utah.edu/common/home/gompert-group2/data/Lycaeides_poolSeq/. The current list of samples follows (this is from a single round of sequencing):

| Population | Samples (N)| Nominal taxon |
|------------|------------|---------------|
| ABM | ABM20 (48) | L. melissa |
| BCR | BCR17 (48), BCR17rep (48) | JH (admixed) |
| BHP | BHP19 (48) | L. melissa |
| BKM | BKM19 (33) | Warners (admixed) |
| BTB | BTB17 (48), BTB17rep (48) | JH (admixed) | 
| CLH | CLH19 (36) | White Mt (admixed) |
| CP | CP19 (48), CP19rep (48) | Sierra (admixed) |
| EP | EP19 (48), EP19rep (48) | Warners (admixed) |
| GNP | GNP17 (56) | L. iads |
| HJ | HJ20 (48) | L. melissa |
| HNV | HNV17 (48) | JH (admixed) |
| LS | LS19 (48) | L. anna |
| MEN | MEN12 (10) | L. argus / argyrognomon (France) |
| MR |  MR20 (48) | Sierra (admixed) |
| MTU | MTU20 (48) | L. melissa |
| SBW | SBW18 (20) | L. idas (Alaska) | 
| SHC | SHC11 (46) | L. anna |
| SIN | SIN10 (48) | L. melissa |
| SUV | SUV20 (51) | L. melissa |
| TBY | TBY11 (24) | L. idas (disjunct population) | 
| TIC | TIC19 (48) | Sierra (admixed?) |
| VE | VE20 (48) | L. melissa |
| YG | YG20 (48) | L. anna |

Replicates (rep) are the same DNA extractions re-pooled and sequenced as a distinct sample. MEN = outgroup. The data were generated and cleaned up by BGI. Reports for each round of sequencing are attached: [BGI_F22FTSUSAT0310-01_LYCgpswR_report_en.pdf](https://github.com/zgompert/LycSpaceTimePoolSeq/files/9940314/BGI_F22FTSUSAT0310-01_LYCgpswR_report_en.pdf)
. 

# DNA Sequence Alignment

I am aligning the DNA sequence data to the updated (based on PacBio) *L. melissa* genome. I am using `bwa-mem2` for this, which is basically just a sped up version of `bwa mem` that also works directly with gzipped files [https://github.com/bwa-mem2/bwa-mem2](https://github.com/bwa-mem2/bwa-mem2). I am using `bwa-mem2` version 2.0pre2. 

First, I (re)indexed the reference genome.

```bash
## index genome with bwa-mem2
/uufs/chpc.utah.edu/common/home/u6000989/source/bwa-mem2-2.0pre2_x64-linux/bwa-mem2 index /uufs/chpc.utah.edu/common/home/gompert-group3/data/LmelGenome/Lmel_dovetailPacBio_genome.fasta
```
Then, I set up the alignment. The submission scrip is (from /uufs/chpc.utah.edu/common/home/gompert-group2/data/Lycaeides_poolSeq/Scripts):

```bash
#!/bin/sh
#SBATCH --time=240:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --account=gompert-kp
#SBATCH --partition=gompert-kp
#SBATCH --job-name=bwa-mem2
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=zach.gompert@usu.edu

module load samtools
##Version: 1.16 (using htslib 1.16)

cd /uufs/chpc.utah.edu/common/home/gompert-group2/data/Lycaeides_poolSeq/Alignments

perl BwaMemFork.pl ../F22FTSUSAT0310-01_LYCgpswR/soapnuke/clean/*/*1.fq.gz 
```

Which runs the following:

```perl
#!/usr/bin/perl
#
# alignment with bwa mem 
#


use Parallel::ForkManager;
my $max = 40;
my $pm = Parallel::ForkManager->new($max);
my $genome = "/uufs/chpc.utah.edu/common/home/gompert-group3/data/LmelGenome/Lmel_dovetailPacBio_genome.fasta";

FILES:
foreach $fq1 (@ARGV){
	$pm->start and next FILES; ## fork
	$fq2 = $fq1;
	$fq2 =~ s/_1\.fq\.gz/_2.fq.gz/ or die "failed substitution for $fq1\n";
        $fq1 =~ m/clean\/([A-Za-z0-9]+)/ or die "failed to match id $fq1\n";
	$ind = $1;
	$fq1 =~ m/([A-Za-z_\-0-9]+)_1\.fq\.gz$/ or die "failed match for file $fq1\n";
	$file = $1;
        system "/uufs/chpc.utah.edu/common/home/u6000989/source/bwa-mem2-2.0pre2_x64-linux/bwa-mem2 mem -t 1 -k 19 -r 1.5 -R \'\@RG\\tID:Lyc-"."$ind\\tLB:Lyc-"."$ind\\tSM:Lyc-"."$ind"."\' $genome $fq1 $fq2 | samtools sort -@ 2 -O BAM -o $ind"."_$file.bam - && samtools index -@ 2 $ind"."_$file.bam\n";

	$pm->finish;
}
```
I am pipping the results on to `samtools` (version 1.16) to compress, sort and index the alignments.

I then merged the bam files for each population using `samtools` version 1.16. This was donw with the following shell script:

```bash
#!/bin/sh
#SBATCH --time=240:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --account=gompert-np
#SBATCH --partition=gompert-np
#SBATCH --job-name=merge
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=zach.gompert@usu.edu

module load samtools
##Version: 1.16 (using htslib 1.16)

cd /uufs/chpc.utah.edu/common/home/gompert-group2/data/Lycaeides_poolSeq/Alignment

perl ../Scripts/MergeFork.pl 
```

Which runs


```perl
#!/usr/bin/perl
#
# merge alignments for each population sample with samtools version XX 
#


use Parallel::ForkManager;
my $max = 40;
my $pm = Parallel::ForkManager->new($max);

open(IDS,"pids.txt");
while(<IDS>){
	chomp;
	push(@IDs,$_);
}
close(IDS);

FILES:
foreach $id (@IDs){
	$pm->start and next FILES; ## fork
        system "samtools merge -c -p -o Merged/$id.bam $id"."_*.bam\n";
	system "samtools index -@ 2 Merged/$id.bam\n";
	$pm->finish;
}

$pm->wait_all_children;
```

pids.txt lists all of the population IDs.

# Removing PCR duplicates

I am using `samtools` (version 1.16) to remove PCR duplicates. This is more efficient than `PicardTools` and performs similarly, see [Ebbert 2016](https://link.springer.com/article/10.1186/s12859-016-1097-3). I am following the standard prtocol from [`samtools`](https://www.htslib.org/doc/samtools-markdup.html). I am using the default option (same as `-m t`) to measure positions based on template start/end. And I am using `-r` to not just mark but remove duplicates. The submission script is:

```bash
#!/bin/sh
#SBATCH --time=240:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --account=gompert-np
#SBATCH --partition=gompert-np
#SBATCH --job-name=dedup
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=zach.gompert@usu.edu

module load samtools
##Version: 1.16 (using htslib 1.16)


cd /scratch/general/nfs1/dedup

perl /uufs/chpc.utah.edu/common/home/gompert-group2/data/Lycaeides_poolSeq/Scripts/RemoveDupsFork.pl *bam
```

Which runs

```perl
#!/usr/bin/perl
#
# PCR duplicate removal with samtools
#


use Parallel::ForkManager;
my $max = 40;
my $pm = Parallel::ForkManager->new($max);

FILES:
foreach $bam (@ARGV){
	$pm->start and next FILES; ## fork
	$bam =~ m/^([A-Za-z0-9]+)/ or die "failed to match $bam\n";
	$base = $1;
	system "samtools collate -o co_$base.bam $bam /scratch/general/nfs1/dedup/t$bam\n";
	system "samtools fixmate -m co_$base.bam fix_$base.bam\n";
	system "samtools sort -o sort_$base.bam fix_$base.bam\n";
	## using default definition of dups
	## measure positions based on template start/end (default). = -m t
	system "markdup -T /scratch/general/nfs1/dedup -r sort_$base.bam dedup_$base.bam\n";
	$pm->finish;
}

$pm->wait_all_children;
```

# Variant Calling

I called variants `bcftools` (version 1.16). I did not perform INDEL realignment (it doesn't seem to be important ) and this approach is not aware I have pooled data. Nonetheless, this should be a solid approach for estimating allele frequencies and quantifying genetic differentiation among pouplations and species. 

For this, I ran the following submission script,

```bash
#!/bin/sh
#SBATCH --time=240:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=24
#SBATCH --account=gompert-np
#SBATCH --partition=gompert-np
#SBATCH --job-name=bcf_call
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=zach.gompert@usu.edu


module load samtools
## version 1.16
module load bcftools
## version 1.16

cd /scratch/general/nfs1/dedup

perl /uufs/chpc.utah.edu/common/home/gompert-group2/data/Lycaeides_poolSeq/Scripts/BcfForkLg.pl chrom*list 
```

which runs

```perl
#!/usr/bin/perl
#
# samtools/bcftools variant calling by LG 
#


use Parallel::ForkManager;
my $max = 26;
my $pm = Parallel::ForkManager->new($max);

my $genome ="/uufs/chpc.utah.edu/common/home/gompert-group3/data/LmelGenome/Lmel_dovetailPacBio_genome.fasta";




foreach $chrom (@ARGV){
	$pm->start and next; ## fork
        $chrom =~ /chrom([0-9\.]+)/ or die "failed here: $chrom\n";
	$out = "o_lycpool_chrom$1";
	system "bcftools mpileup -b bams -d 1000 -f $genome -R $chrom -a FORMAT/DP,FORMAT/AD -q 20 -Q 30 -I -Ou | bcftools call -v -c -p 0.01 -Ov -o $out"."vcf\n";
	#system "bcftools mpileup -S bams -d 1000 -f $genome -R $chrom -q 20 -Q 30 -I -t DP,AD,ADF,ADR -Ou | bcftools call -v -c -p 0.01 -Ov -o $out.vcf\n";


	$pm->finish;

}

$pm->wait_all_children;
```
Note that each chromosome (big scaffold) is being processes separately (chrom*list). 

# Preliminary analyses based on `bcftools` results

I did some quick, preliminary analyses based on the above calling with `bcftools`. All of this is in /uufs/chpc.utah.edu/common/home/gompert-group2/data/Lycaeides_poolSeq/Variants/.

First, I used `GATK` version (4.1.4.1) for some simple filtering of the variants, keeping only those with mapping quality > 30 and a depth > 1750.

```perl
#!/usr/bin/perl
#
# filter vcf files 
#


use Parallel::ForkManager;
my $max = 48;
my $pm = Parallel::ForkManager->new($max);


foreach $vcf (@ARGV){
	$pm->start and next; ## fork
	$o = $vcf;
	$o =~ s/o_// or die "failed sub $vcf\n";
	system "bgzip $vcf\n";
	system "tabix $vcf.gz\n";
	system "java -jar /uufs/chpc.utah.edu/sys/installdir/gatk/gatk-4.1.4.1/gatk-package-4.1.4.1-local.jar IndexFeatureFile -I $vcf.gz\n";
	system "java -jar /uufs/chpc.utah.edu/sys/installdir/gatk/gatk-4.1.4.1/gatk-package-4.1.4.1-local.jar VariantFiltration -R /uufs/chpc.utah.edu/common/home/gompert-group3/data/LmelGenome/Lmel_dovetailPacBio_genome.fasta -V $vcf.gz -O filt_$o.gz --filter-name \"depth\" --filter-expression \"DP < 1750\" --filter-name \"mapping\" --filter-expression \"MQ < 30\"\n";
	system "bgzip -d filt_$o.gz\n";

	$pm->finish;

}

$pm->wait_all_children;
```

Then I extracted the allele depths from the filtered vcf files,

```bash
#!/usr/bin/bash
#
# extract allele depth AD from biallelic SNPs that passed filtering 
#

for f in filt*vcf
do
	echo "Processing $f"
	out="$(echo $f | sed -e 's/vcf/txt/')"
	echo "Output is ad1_$out"
	grep ^Sc $f | grep PASS | grep -v [ATCG],[ATCG] | perl -p -i -e 's/^.+AD\s+//' | perl -p -i -e 's/\S+:(\d+),(\d+)/\1/g' > ad1_$out   
	grep ^Sc $f | grep PASS | grep -v [ATCG],[ATCG] | perl -p -i -e 's/^.+AD\s+//' | perl -p -i -e 's/\S+:(\d+),(\d+)/\2/g' > ad2_$out
done
```
This creates allele depth files for each allele (ad1* and ad2*) and chromosome, which I am using to obtain maximum likelihood allele frequency estimates (this is just the sample frequency... even a beta-binomial model would give the same point estimate with a flat prior). I conducted PCAs and estimates of Fst by chromosome, along with window-based Fst scans for some populations. The code is [pcaFst.R](pcaFst.R). It all generally makes sense. The Eurasian sample is clearly the outgroup (followed by the Alaskan samples), replicates are very similar, and so are temporal samples to a lesser extent. Most chromosomes have similar samples though the Z stands out, and interstingly puts ABM closer to idas and anna. The PCA is attached [WG_LG_PCAs.pdf](https://github.com/zgompert/LycSpaceTimePoolSeq/files/10094994/WG_LG_PCAs.pdf).


