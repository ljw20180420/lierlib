from lierlib import *

### mmej probability predict
for M in range(1,6):
    longest_common_substring_cumulative(f"/home/ljw/wuqiang/test_run/{M}.mmej.table", n=numpy.arange(9,101), M=M, p=0.25)

### read map analysis
genome_file = '/home/ljw/hg19_with_bowtie2_index/hg19.fa'
excu = '/home/ljw/new_fold/old_desktop/shoujia/Graph_Projects/GAP/build/GAP'
fastq_files = ['/home/ljw/wuqiang/test_run/high_quality_cutadapt_barcode_WTdw2-STM_R1.fq']
chrss1 = [['chr11','chr11','chr11']]
chrss2 = [['chr11','chr11','chr11']]
cutss1 = [[4093764,4093764,4093764]]
cutss2 = [[4093764,4114684,4114684]]
strandss1 = [['+','+','+']]
strandss2 = [['+','+','-']]
threads = 6

best_files, count_files = [], []
for fastq_file, chrs1, chrs2, cuts1, cuts2, strands1, strands2 in zip(fastq_files, chrss1, chrss2, cutss1, cutss2, strandss1, strandss2):
    count_file, fasta_file = count_fastq(fastq_file)
    count_files.append(count_file)

    # this step is slow
    align(fasta_file, genome_file, excu, msn=3, mr=4, T=-10, dT=-5, bs=10, tz=threads)

    to_sam(fasta_file)
    subprocess.check_output(f'samtools view -@ {threads} -b -T {genome_file} {f"{fasta_file}.sam"} | samtools sort -@ {threads} -o {f"{fasta_file}.bam"}', shell=True)
    subprocess.check_output(f'samtools index -b {f"{fasta_file}.bam"}', shell=True)

    # # not used, do not excute
    # blocks = []
    # for chr1, chr2, cut1, cut2, strand1, strand2 in zip(chrs1, chrs2, cuts1, cuts2, strands1, strands2):
    #     if chr1!=chr2:
    #         raise Exception("chrom1 and chrom2 must be the same")
    #     blocks.append(filter_bam(f"{fasta_file}.bam", (chr1, cut1-100000, cut2+100000), [strand1,strand2], ord='forward', gap=100, segnums=[2]))

    indel_files = []
    for chr1, chr2, cut1, cut2, strand1, strand2 in zip(chrs1, chrs2, cuts1, cuts2, strands1, strands2):
        indel_files.append(get_indel(fasta_file, count_file, genome_file, (chr1, cut1-110000, cut1+110000), strand1, cut1, (chr2, cut2-110000, cut2+110000), strand2, cut2, thres=100000, align_mode="up"))
    best_files.append(f'{fasta_file}.indel')
    best_indel(indel_files, best_files[-1])

    dis_range = [min(cuts1+cuts2)-1000, max(cuts1+cuts2)+10000]
    draw_coverage(best_files[-1], chrs1, cuts1, strands1, chrs2, cuts2, strands2, dis_range, align_mode='align_left')

    cut12 = list(set(cuts1 + cuts2))
    if len(cut12)!=2:
        raise Exception("total cut positions must be two")
    for chr1, chr2, cut1, cut2, strand1, strand2 in zip(chrs1, chrs2, cuts1, cuts2, strands1, strands2):
        for second_resection_range in [[0,99],[100,numpy.inf]]:
            draw_large_right_resection(best_files[-1], chr1, cut1, strand1, chr2, cut2, strand2, dis_range, cut12, second_resection_range=second_resection_range, align_mode='align_left')
        draw_segs(best_files[-1], chr1, cut1, strand1, chr2, cut2, strand2, dis_range, cut12, align_mode='align_left')
        count_types(best_files[-1], count_file, chr1, cut1, strand1, chr2, cut2, strand2)
    
    if all([chr1==chr2 for chr1, chr2 in zip(chrs1, chrs2)]):
        subprocess.check_output(f'samtools tview -p {f"{chrs1[0]}:{dis_range[0]}"} -d T -w {dis_range[1]-dis_range[0]} {f"{fasta_file}.bam"} {genome_file} > {fasta_file}.display', shell=True)

### this runs out of loop
summary_count(best_files, count_files, chrss1, chrss2, cutss1, cutss2, strandss1, strandss2)

    