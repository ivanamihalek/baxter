import os

from bad_bac_exercise.scripts.utils import die_if_not_executable, die_if_not_nonempty_file, run_subprocess, randomword


def bwa_mem2_alignment(reference_genome: str, read_files: list[str], out_bam: str):
    # todo - check the sanity of out_bam
    # is the directory writable? does the name start with "." ?

    # first need to index the reference
    # todo - this should be in some settings file or some such
    aligner  = "/usr/third/bwa-mem2-2.2.1/bwa-mem2"
    samtools = "/usr/third/samtools-1.15.1/samtools"
    for executable in [aligner, samtools]:
        die_if_not_executable(executable)

    for input_file in [reference_genome, *read_files]:
        die_if_not_nonempty_file(input_file)

    # create index file, if not present
    if not os.path.exists(f"{reference_genome}.bwt.2bit.64"):
        # /usr/third/bwa-mem2-2.2.1/bwa-mem2 index GCF_000195955.2.fa
        cmd = f"{aligner} index {reference_genome}"
        ret = run_subprocess(cmd)

    # align
    # /usr/third/bwa-mem2-2.2.1/bwa-mem2 mem GCF_000195955.2.fa sample_reads_R1.fastq sample_reads_R2.fastq  > test.sam
    cmd = " ".join([f"{aligner} mem {reference_genome}"] + read_files)
    tmp_sam_name = f"{randomword(5)}.{os.getpid()}.sam"
    ret = run_subprocess(cmd, stdoutfnm=tmp_sam_name)

    # sort the output in sam format, and turn it into bam
    # /usr/third/samtools-1.15.1/samtools sort test.sam > test.sorted.bam
    bam_name = out_bam if out_bam[:-4] == ".bam" else f"{out_bam}.bam"
    cmd = f"{samtools} sort {tmp_sam_name}"
    ret = run_subprocess(cmd, stdoutfnm=bam_name)
    # remove the original sam file
    os.remove(tmp_sam_name)

    # create bam index
    # /usr/third/samtools-1.15.1/samtools index test.sorted.bam
    cmd = f"{samtools} index {bam_name}"
    ret = run_subprocess(cmd)


def gatk_haplotyper_variant_caller():
    # implement variant calling to check where didi the variants end up
    # https://gatk.broadinstitute.org/hc/en-us/articles/360037225632-HaplotypeCaller
    # java -jar /usr/third/gatk/picard.jar AddOrReplaceReadGroups I=test.sorted.bam O=test.sorted.read_groups.bam SORT_ORDER=coordinate RGID=foo RGLB=bar  RGPL=illumina RGSM=Sample1 CREATE_INDEX=True RGPU=unit1
    # /usr/third/gatk/gatk-4.6.1.0/gatk CreateSequenceDictionary -R GCF_000195955.2.fa
    # fai also must be present
    # /usr/third/gatk/gatk-4.6.1.0/gatk --java-options "-Xmx4g" HaplotypeCaller -R GCF_000195955.2.fa  -I test.sorted.read_groups.bam -O tes.out.vcf
    #
    # cmd = "/usr/third/bwa-mem2-2.2.1/bwa-mem2 "
    # ret = run_subprocess(cmd)
    pass
