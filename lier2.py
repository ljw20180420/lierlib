from lierlib import *

### read map analysis
genome_file = '/home/ljw/hg19_with_bowtie2_index/hg19.fa'
excu = '/home/ljw/new_fold/old_desktop/shoujia/projects/Rearrangement/build/rearrangement'
fastq_files = ['/home/ljw/wuqiang/test_run2/high_quality_cutadapt_barcode_WTdw2-STM_R1.fq']
chrss1 = [['chr11','chr11','chr11']]
chrss2 = [['chr11','chr11','chr11']]
cutss1 = [[4093764,4093764,4093764]]
cutss2 = [[4093764,4114684,4114684]]
strandss1 = [['+','+','+']]
strandss2 = [['+','+','-']]
threads = 12
refext = 11000
inner, outter = 100000, 110000
minseg = 20
absTini = 10

best_files, count_files = [], []
for fastq_file, chrs1, chrs2, cuts1, cuts2, strands1, strands2 in zip(fastq_files, chrss1, chrss2, cutss1, cutss2, strandss1, strandss2):
    # simple check
    if not (chrs1[0]==chrs1[1] and chrs1[1]==chrs1[2]):
        raise Exception("chrs1 must all be the same")
    cut12 = list(set(cuts1 + cuts2))
    if len(cut12)!=2:
        raise Exception("total cut positions must be two")

    # unify reads
    count_file, fasta_file = count_fastq(fastq_file)
    count_files.append(count_file)

    # align all possible references
    for chr1, chr2, cut1, cut2, strand1, strand2, i in zip(chrs1+[chrs1[0]], chrs2+[None], cuts1+[cuts1[0]], cuts2+[None], strands1+[strands1[0]], strands2+[None], numpy.arange(len(chrs1)+1)):
        shutil.copy(count_file, f"{count_file}.{i}")
        # token = 1700
        # with open(f"{count_file}.{i}", "w") as fw, open(count_file, "r") as fr:
        #     for rn in range(token):
        #         _ = fw.write(fr.readline())
        rearrangement(f"{count_file}.{i}", threads, excu, genome_file, refext, chr1=chr1, cut1=cut1, strand1=strand1, chr2=chr2, cut2=cut2, strand2=strand2)
        with open(f"{count_file}.{i}.alg", 'w') as fw:
            j = 0
            while os.path.exists(f"{os.path.join(os.path.dirname(count_file),'tmp',os.path.basename(count_file))}.{i}.{j}.alg"):
                with open(f"{os.path.join(os.path.dirname(count_file),'tmp',os.path.basename(count_file))}.{i}.{j}.alg", 'r') as fr:
                    fw.writelines(fr.readlines())
                j += 1
    shutil.rmtree(os.path.join(os.path.dirname(count_file),'tmp'))    

    # select the best reference
    with open(f"{count_file}.alg", "w") as fw, open(f"{count_file}.0.alg", "r") as fr0, open(f"{count_file}.1.alg", "r") as fr1, open(f"{count_file}.2.alg", "r") as fr2, open(f"{count_file}.3.alg", "r") as fr3:
        jbs = []
        while True:
            liness = [[fr.readline() for i in range(3)] for fr in [fr0, fr1, fr2, fr3]]
            if not liness[0][0]:
                break
            scores = numpy.array([int(lines[0].split('\t')[2]) for lines in liness])
            # for the single cut case, we give it a bounus of absTini
            scores[-1] += absTini
            jb = numpy.argmax(scores) 
            _ = fw.writelines(liness[jb])
            jbs.append(jb)
    rearrange_to_GAP(fasta_file, f"{count_file}.alg", jbs, chrs1+[chrs1[0]], cuts1+[cuts1[0]], strands1+[strands1[0]], chrs2+[None], cuts2+[None], strands2+[None], genome_file, refext, minseg=minseg)

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
        indel_files.append(get_indel(fasta_file, count_file, genome_file, (chr1, cut1-outter, cut1+outter), strand1, cut1, (chr2, cut2-outter, cut2+outter), strand2, cut2, thres=inner, align_mode="up"))
    best_files.append(f'{fasta_file}.indel')
    best_indel(indel_files, best_files[-1])

    dis_range = [min(cuts1+cuts2)-1000, max(cuts1+cuts2)+10000]
    draw_coverage(best_files[-1], chrs1, cuts1, strands1, chrs2, cuts2, strands2, dis_range, align_mode='align_left')

    for chr1, chr2, cut1, cut2, strand1, strand2 in zip(chrs1, chrs2, cuts1, cuts2, strands1, strands2):
        # for second_resection_range in [[0,99],[100,numpy.inf]]:
        #     draw_large_right_resection(best_files[-1], chr1, cut1, strand1, chr2, cut2, strand2, dis_range, cut12, second_resection_range=second_resection_range, align_mode='align_left')
        # draw_segs(best_files[-1], chr1, cut1, strand1, chr2, cut2, strand2, dis_range, cut12, align_mode='align_left')
        count_types(best_files[-1], count_file, chr1, cut1, strand1, chr2, cut2, strand2)
    
    if all([chr1==chr2 for chr1, chr2 in zip(chrs1, chrs2)]):
        subprocess.check_output(f'samtools tview -p {f"{chrs1[0]}:{dis_range[0]}"} -d T -w {dis_range[1]-dis_range[0]} {f"{fasta_file}.bam"} {genome_file} > {fasta_file}.display', shell=True)

### this runs out of loop
summary_count(best_files, count_files, chrss1, chrss2, cutss1, cutss2, strandss1, strandss2)

### mmej probability predict
for M in range(1,6):
    longest_common_substring_cumulative(f"/home/ljw/wuqiang/test_run/{M}.mmej.table", n=numpy.arange(9,101), M=M, p=0.25)

