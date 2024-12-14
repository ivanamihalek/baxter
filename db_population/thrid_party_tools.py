import os

from db_population.utils import die_if_not_executable, die_if_not_nonempty_file, run_subprocess, randomword


def bwa_mem2_alignment(reference_genome: str, workdir: str, read_files: list[str], out_bam: str):
    # todo - check the sanity of out_bam
    # is the directory writable? does the name start with "." ?
    i_came_from = os.getcwd()
    os.chdir(workdir)

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
    bam_name = out_bam if out_bam[-4:] == ".bam" else f"{out_bam}.bam"
    cmd = f"{samtools} sort {tmp_sam_name} -o bam_name"
    ret = run_subprocess(cmd)
    die_if_not_nonempty_file(bam_name)
    # remove the original sam file
    os.remove(tmp_sam_name)

    # create bam index
    # /usr/third/samtools-1.15.1/samtools index test.sorted.bam
    cmd = f"{samtools} index {bam_name}"
    ret = run_subprocess(cmd)

    os.chdir(i_came_from)




def gatk_haplotyper_variant_caller(reference_genome: str, workdir: str,  in_bam: str, out_vcf: str):
    i_came_from = os.getcwd()
    os.chdir(workdir)

    # implement variant calling to check where did the variants end up
    # https://gatk.broadinstitute.org/hc/en-us/articles/360037225632-HaplotypeCaller
    java = "/usr/lib/jvm/java-21-openjdk-amd64/bin/java"
    gatk = "/usr/third/gatk/gatk-4.6.1.0/gatk"
    samtools = "/usr/third/samtools-1.15.1/samtools"
    for executable in [java, gatk, samtools]:
        die_if_not_executable(executable)

    picard_jar = "/usr/third/gatk/picard.jar"
    for input_file in [picard_jar, reference_genome, in_bam]:
        die_if_not_nonempty_file(input_file)

    if in_bam[-4:] != ".bam":
        print(f"the variant calling input file {in_bam} does not seem to have the .bam extension")
        exit()

    # java -jar /usr/third/gatk/picard.jar AddOrReplaceReadGroups I=test.sorted.bam O=test.sorted.read_groups.bam
    # SORT_ORDER=coordinate RGID=foo RGLB=bar  RGPL=illumina RGSM=Sample1 CREATE_INDEX=True RGPU=unit1
    read_grps_bam = in_bam[:-4] + ".read_groups.bam"
    cmd = f"{java} -jar {picard_jar} AddOrReplaceReadGroups I={in_bam} O={read_grps_bam} "
    cmd += "SORT_ORDER=coordinate RGID=foo RGLB=bar  RGPL=illumina RGSM=Sample1 CREATE_INDEX=True RGPU=unit1"
    ret = run_subprocess(cmd)

    if reference_genome[-3:] != ".fa":
        print(f"the reference genome {reference_genome} does not seem to have the .fa extension")
        exit()

    dict_file = reference_genome[:-3] + ".dict"
    if not os.path.exists(dict_file):
        # /usr/third/gatk/gatk-4.6.1.0/gatk CreateSequenceDictionary -R GCF_000195955.2.fa
        cmd = f"{gatk}  CreateSequenceDictionary -R {reference_genome}"
        ret = run_subprocess(cmd)

    # fai also must be present
    if not os.path.exists(reference_genome + ".fai"):
        cmd = f"{samtools}  faidx  {reference_genome}"
        ret = run_subprocess(cmd)

    # and finally, call the variants
    # /usr/third/gatk/gatk-4.6.1.0/gatk --java-options "-Xmx4g" HaplotypeCaller -R GCF_000195955.2.fa
    # -I test.sorted.read_groups.bam -O tes.out.vcf
    cmd = f"{gatk} HaplotypeCaller -R {reference_genome} -I {read_grps_bam} -O {out_vcf}"
    ret = run_subprocess(cmd)

    die_if_not_nonempty_file(out_vcf, prepend="after gatk variant calling step")
    os.chdir(i_came_from)
