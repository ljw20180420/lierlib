from lierlib import *
import requests

### download 
# file='https://hgdownload.cse.ucsc.edu/goldenpath/hg19/database/rmsk.txt.gz'
# r = requests.get(file)
# with open(os.path.basename(file), 'wb') as f:
#     f.write(r.content)

### read map analysis
genome_file = '/home/ljw/hg19_with_bowtie2_index/hg19.fa'
excu = '/home/ljw/wuqiang/lierlib/Rearrangement/build/rearrangement'
datadir = '/home/ljw/wuqiang/sxdata'
fastq_files_dict = {}
fastq_files = []
for sample in ['DCK','HPRT']:
    for primer in ['U','D']:
        for Ri in ['R1', 'R2']:
            fastq_files_dict[f'{sample}_{primer}_{Ri}'] = [os.path.join(datadir,file) for file in os.listdir(datadir) if os.path.isfile(os.path.join(datadir,file)) and file.find(f'-{sample}-{primer}-')!=-1 and file.endswith(f".{Ri}.fastq.gz")]
            fastq_files += fastq_files_dict[f'{sample}_{primer}_{Ri}']

rmsk = bioframe.read_table('rmsk.txt.gz')[[5, 6, 7]]
rmsk.columns = ['chrom', 'start', 'end']

cts = [{'sample' : 'HPRT', 'chr' : 'chrX', 'sgrnaU' : "GTAAGTAAGATCTTAAAATGAGG", 'sgrnaD' : "CCTGATTTTATTTCTGTAGGACT"}, {'sample' : 'DCK', 'chr' : 'chr4', 'sgrnaU' : "ATATTTAGAACTCTTTTCAGTGG", 'sgrnaD' : "CCCTGCCTTTTTCTTCCATCTCT"}]

genome = bioframe.load_fasta(genome_file)
for ct in cts:
    chrgenome = genome[ct['chr']].ff.fetch(ct['chr'])
    ct['cutU'] = chrgenome.find(ct['sgrnaU'])+len(ct['sgrnaU'])-6
    ct['cutD'] = chrgenome.find(ct['sgrnaD'])+6
    ct['rmsk'] = bioframe.select(rmsk, (ct['chr'], ct['cutU']-3000, ct['cutD']+3000))
    # ct['crmsk'] = pandas.DataFrame({'chrom' : [ct['chr']] * (len(ct['rmsk'])+1), 'start' : [ct['cutU']-3000] + ct['rmsk']['end'].values.tolist(), 'end' : ct['rmsk']['start'].values.tolist() + [ct['cutD']+3000]})

threads = 24
inner, outter = 100000, 110000
minseg = 20
absTini = 10
ext_file = [os.path.join(datadir,file) for file in os.listdir(datadir) if os.path.isfile(os.path.join(datadir,file)) and file.endswith(f".fastq.gz.fa.ext")]
for fastq_file in fastq_files:
    if fastq_file.endswith(".R2.fastq.gz"):
        continue
    
    # if f"{fastq_file}.fa.ext" in ext_file:
    #     continue

    for ct in cts:
        if fastq_file.find(ct['sample'])!=-1:
            sample = ct['sample']
            chrs1, chrs2, cutU, cutD = [ct['chr']]*3, [ct['chr']]*3, ct['cutU'], ct['cutD']
            refext = ct['cutD'] - ct['cutU'] + 150
            ctt = ct
            break
    
    if fastq_file.find(f'-{sample}-U-')!=-1:
        cuts1, cuts2, strands1, strands2 = [cutU]*3, [cutU,cutD,cutD], ['+','+','+'], ['+','+','-']
    else:
        cuts1, cuts2, strands1, strands2 = [cutD]*3, [cutD,cutU,cutU], ['-','-','-'], ['-','-','+']

    # simple check
    if not (chrs1[0]==chrs1[1] and chrs1[1]==chrs1[2]):
        raise Exception("chrs1 must all be the same")
    cut12 = list(set(cuts1 + cuts2))
    if len(cut12)!=2:
        raise Exception("total cut positions must be two")

    # unify reads
    # count_file, fasta_file = count_fastq(fastq_file)
    count_file, fasta_file = f'{fastq_file}.txt', f'{fastq_file}.fa'
    best_file = f'{fasta_file}.indel'

    # # align all possible references
    # for chr1, chr2, cut1, cut2, strand1, strand2, i in zip(chrs1+[chrs1[0]], chrs2+[None], cuts1+[cuts1[0]], cuts2+[None], strands1+[strands1[0]], strands2+[None], numpy.arange(len(chrs1)+1)):
    #     shutil.copy(count_file, f"{count_file}.{i}")
    #     rearrangement(f"{count_file}.{i}", threads, excu, genome_file, refext, chr1=chr1, cut1=cut1, strand1=strand1, chr2=chr2, cut2=cut2, strand2=strand2)
    #     with open(f"{count_file}.{i}.alg", 'w') as fw:
    #         j = 0
    #         while os.path.exists(f"{os.path.join(os.path.dirname(count_file),'tmp',os.path.basename(count_file))}.{i}.{j}.alg"):
    #             with open(f"{os.path.join(os.path.dirname(count_file),'tmp',os.path.basename(count_file))}.{i}.{j}.alg", 'r') as fr:
    #                 fw.writelines(fr.readlines())
    #             j += 1
    #     os.remove(f"{count_file}.{i}")
    # shutil.rmtree(os.path.join(os.path.dirname(count_file),'tmp'))    

    # # select the best reference
    # with open(f"{count_file}.alg", "w") as fw, open(f"{count_file}.0.alg", "r") as fr0, open(f"{count_file}.1.alg", "r") as fr1, open(f"{count_file}.2.alg", "r") as fr2, open(f"{count_file}.3.alg", "r") as fr3:
    #     jbs = []
    #     while True:
    #         liness = [[fr.readline() for i in range(3)] for fr in [fr0, fr1, fr2, fr3]]
    #         if not liness[0][0]:
    #             break
    #         scores = numpy.array([int(lines[0].split('\t')[2]) for lines in liness])
    #         # for the single cut case, we give it a bounus of absTini
    #         scores[-1] += absTini
    #         jb = numpy.argmax(scores) 
    #         _ = fw.writelines(liness[jb])
    #         jbs.append(jb)
    # for i in range(4):
    #     os.remove(f"{count_file}.{i}.alg")
    # rearrange_to_GAP(fasta_file, f"{count_file}.alg", jbs, chrs1+[chrs1[0]], cuts1+[cuts1[0]], strands1+[strands1[0]], chrs2+[None], cuts2+[None], strands2+[None], genome_file, refext, minseg=minseg)
    # os.remove(f"{count_file}.alg")

    # to_sam(fasta_file)
    # subprocess.check_output(f'samtools view -@ {threads} -b -T {genome_file} {f"{fasta_file}.sam"} | samtools sort -@ {threads} -o {f"{fasta_file}.bam"}', shell=True)
    # subprocess.check_output(f'samtools index -b {f"{fasta_file}.bam"}', shell=True)

    # # # not used, do not excute
    # # blocks = []
    # # for chr1, chr2, cut1, cut2, strand1, strand2 in zip(chrs1, chrs2, cuts1, cuts2, strands1, strands2):
    # #     if chr1!=chr2:
    # #         raise Exception("chrom1 and chrom2 must be the same")
    # #     blocks.append(filter_bam(f"{fasta_file}.bam", (chr1, cut1-100000, cut2+100000), [strand1,strand2], ord='forward', gap=100, segnums=[2]))

    # indel_files = []
    # for chr1, chr2, cut1, cut2, strand1, strand2 in zip(chrs1, chrs2, cuts1, cuts2, strands1, strands2):
    #     indel_files.append(get_indel(fasta_file, count_file, genome_file, (chr1, cut1-outter, cut1+outter), strand1, cut1, (chr2, cut2-outter, cut2+outter), strand2, cut2, thres=inner, align_mode="up"))
    # best_indel(indel_files, best_file)

    dis_range = [min(cuts1+cuts2)-refext-300, max(cuts1+cuts2)+refext+300]
    rmsk_thres = 1 # segment2 with rmsk ratio larger than rmsk_thres will be filtered
    # subprocess.check_output(f"mv {best_file} {best_file}.bak", shell=True)
    tmp_best_indel = pandas.read_csv(f"{best_file}.bak", sep='\t')
    tmp_best_indel_ts = tmp_best_indel[tmp_best_indel['ref2']!='*'][['ref2','mapstart2','mapend2']].rename(columns={'ref2' : 'chrom', 'mapstart2' : 'start', 'mapend2' : 'end'})
    tmp_best_indel_ts[['start', 'end']] = tmp_best_indel_ts[['start', 'end']].astype(int)
    ts_index = tmp_best_indel_ts.index
    tmp_best_indel_ts = bioframe.coverage(tmp_best_indel_ts, ctt['rmsk'])
    tmp_best_indel['map2_not_rmsk'] = True
    tmp_best_indel['map2_not_rmsk'][ts_index] = [True if coverage/(end-start)<rmsk_thres else False for start, end, coverage in zip(tmp_best_indel_ts['start'], tmp_best_indel_ts['end'], tmp_best_indel_ts['coverage'])]
    tmp_best_indel['key1_prop'] = (tmp_best_indel['strand1']=="+") & (tmp_best_indel['key1']<=tmp_best_indel['cut1']) & (tmp_best_indel['key1']>=tmp_best_indel['cut1']-130) | (tmp_best_indel['strand1']=="-") & (tmp_best_indel['key1']>=tmp_best_indel['cut1']) & (tmp_best_indel['key1']<=tmp_best_indel['cut1']+130)
    tmp_best_indel['map1_long'] = (tmp_best_indel['mapstart1']=="*") | (tmp_best_indel['mapend1']-tmp_best_indel['mapstart1']>=25)
    tmp_best_indel['map2_long'] = [True if mapstart2=="*" or int(mapend2)-int(mapstart2)>=25 else False for mapstart2, mapend2 in zip(tmp_best_indel['mapstart2'], tmp_best_indel['mapend2'])]
    # tmp_best_indel[tmp_best_indel['map2_not_rmsk'] & tmp_best_indel['key1_prop'] & tmp_best_indel['map1_long'] & tmp_best_indel['map2_long']].reset_index(drop=True).to_csv(f'{best_file}', header=True, index=False, sep='\t')
    tmp_best_indel[tmp_best_indel['map2_not_rmsk'] & tmp_best_indel['key1_prop'] & tmp_best_indel['map2_long']].reset_index(drop=True).to_csv(f'{best_file}', header=True, index=False, sep='\t')

    draw_coverage(best_file, chrs1, cuts1, strands1, chrs2, cuts2, strands2, dis_range, align_mode='align_left')

    for chr1, chr2, cut1, cut2, strand1, strand2 in zip(chrs1, chrs2, cuts1, cuts2, strands1, strands2):
        for second_resection_range in [[0,99],[100,numpy.inf]]:
            draw_large_right_resection(best_file, chr1, cut1, strand1, chr2, cut2, strand2, dis_range, cut12, second_resection_range=second_resection_range, align_mode='align_left')
        draw_segs(best_file, chr1, cut1, strand1, chr2, cut2, strand2, dis_range, cut12, align_mode='align_left')
        count_types(best_file, count_file, chr1, cut1, strand1, chr2, cut2, strand2)
    
    if all([chr1==chr2 for chr1, chr2 in zip(chrs1, chrs2)]):
        subprocess.check_output(f'samtools tview -p {f"{chrs1[0]}:{dis_range[0]}"} -d T -w {dis_range[1]-dis_range[0]} {f"{fasta_file}.bam"} {genome_file} > {fasta_file}.display', shell=True)



### this runs out of loop
best_files, count_files, chrss1, chrss2, cutss1, cutss2, strandss1, strandss2 = [], [], [], [], [], [], [], []
for fastq_file in fastq_files:
    if fastq_file.endswith(".R2.fastq.gz"):
        continue
    
    for ct in cts:
        if fastq_file.find(ct['sample'])!=-1:
            sample = ct['sample']
            chrs1, chrs2, cutU, cutD = [ct['chr']]*3, [ct['chr']]*3, ct['cutU'], ct['cutD']
            refext = ct['cutD'] - ct['cutU'] + 150
            break
    
    if fastq_file.find(f'-{sample}-U-')!=-1:
        cuts1, cuts2, strands1, strands2 = [cutU]*3, [cutU,cutD,cutD], ['+','+','+'], ['+','+','-']
    else:
        cuts1, cuts2, strands1, strands2 = [cutD]*3, [cutD,cutU,cutU], ['-','-','-'], ['-','-','+']

    chrss1.append(chrs1)
    chrss2.append(chrs2)
    cutss1.append(cuts1)
    cutss2.append(cuts2)
    strandss1.append(strands1)
    strandss2.append(strands2)

    count_files.append(f"{fastq_file}.txt")
    best_files.append(f"{fastq_file}.fa.indel")

summary_count(best_files, count_files, chrss1, chrss2, cutss1, cutss2, strandss1, strandss2)

### mmej probability predict
for M in range(1,6):
    longest_common_substring_cumulative(f"/home/ljw/wuqiang/test_run/{M}.mmej.table", n=numpy.arange(9,101), M=M, p=0.25)

### generate randome DNA
import random 
dna_bases = ['A', 'T', 'C', 'G']
dna = ''
n = 100
for i in range(n):
    dna += random.choice(dna_bases)
