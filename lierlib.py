import os, collections, subprocess, pandas, bioframe, itertools, Bio.Seq, pysam, numpy, more_itertools, matplotlib.pyplot, pdb, seaborn, shutil, re, warnings, gzip

### count read numbers and sort by numbers
def count_fastq(fastq_file):
    with gzip.open(fastq_file, 'r') as f:
        duplines = collections.Counter(f.readlines()[1::4])
        with open(f'{fastq_file}.txt', 'w') as fw, open(f'{fastq_file}.fa', 'w') as fa:
            for i, seq_c in enumerate(duplines.most_common(), start=1):
                dec_seq = seq_c[0].decode()
                _ = fw.write(f"{dec_seq[:-1]}\t{seq_c[1]}\n")
                _ = fa.write(f'>seq{i}\n{dec_seq}')
    
    return f'{fastq_file}.txt', f'{fastq_file}.fa'

### use the GAP to align reads
def align(fasta_file, genome_file, excu, msn=1, mr=4, T=-10, dT=-5, bs=10, tz=24):
    run=os.path.basename(fasta_file)
    cwd = os.getcwd()
    os.chdir(os.path.dirname(fasta_file))
    subprocess.check_output(f'{excu} ---threads_sz {tz} ---run {run} ---reads {fasta_file} ---max_round {mr} ---min_seg_num 1 ---max_seg_num {msn} ---block_size {bs} ---nodes --names node0 ---roots node0 ---targets node0 ---globals --names {genome_file} --tails node0 --heads node0 --Ts {T} --dTs {dT} --r true', shell=True)
    if os.path.exists(f'{run}{mr}.fail'):
        subprocess.check_output(f'{excu} ---threads_sz {tz} ---run {run}_last ---reads {run}{mr}.fail ---block_size {bs} ---nodes --names {" ".join([f"node{i}" for i in range(msn+1)])} ---roots node0 ---targets node{msn} ---globals --names {" ".join([genome_file]*msn)} --tails {" ".join([f"node{i}" for i in range(msn)])} --heads {" ".join([f"node{i}" for i in range(1,msn+1)])} --Ts 0 --r {" ".join(["true"]*msn)}', shell=True)
    os.chdir(cwd)

### use the rearrangement to align reads
def rearrangement(count_file, threads, excu, genome_file, refext, chr1=None, cut1=None, strand1=None, chr2=None, cut2=None, strand2=None):
    genome = bioframe.load_fasta(genome_file)

    cwd = os.getcwd()
    os.chdir(os.path.dirname(count_file))
    with open("reference.ref", 'w') as f:
        ref1 = genome[chr1].ff.fetch(chr1, cut1-refext, cut1+refext)
        if strand1=='-':
            ref1 = Bio.Seq.Seq(ref1).reverse_complement().__str__()
        _ = f.write(f"0\t0\t0\n{ref1}\n{refext}\t{refext}\t{refext}\n")
        if chr2!=None:
            ref2 = genome[chr2].ff.fetch(chr2, cut2-refext, cut2+refext)
            if strand2=='-':
                ref2 = Bio.Seq.Seq(ref2).reverse_complement().__str__()
            _ = f.write(f"{refext}\t{refext}\t{refext}\n{ref2}\n{refext*2}\t{refext*2}\t{refext*2}\n")
    subprocess.check_output(f"{excu} -file {count_file} -ref_file reference.ref -ALIGN_MAX 1 -THR_NUM {threads}", shell=True)
    os.chdir(cwd)

### reformat the output of rearrange to that of GAP
def rearrange_to_GAP(fasta_file, best_alg, jbs, chrs1, cuts1, strands1, chrs2, cuts2, strands2, genome_file, refext, minseg=20):
    ref2len = {}
    for row in bioframe.read_table('/home/ljw/hg19_with_bowtie2_index/hg19.chrom.sizes', schema=['chrom', 'length']).itertuples():
        ref2len[row.chrom] = row.length
    # reformat
    with open(best_alg, 'r') as fr, open(f"{fasta_file}.alg", 'w') as fw, open(f"{fasta_file}.ext", 'w') as ext:
        jj=-1
        for info, refam, queryam in more_itertools.batched(fr, 3):
            jj+=1
            info_sp = info.split('\t')
            id, mids, starts, ends = int(info_sp[0]), [info_sp[3],info_sp[6]], [int(info_sp[4])], [int(info_sp[5])]
            if jbs[jj]<len(chrs1)-1:
                starts.append(int(info_sp[7]))
                ends.append(int(info_sp[8]))
            lowidxs = [-1]
            lowtimes = 2 if jbs[jj]==len(chrs1)-1 else 4
            for i in range(lowtimes):
                lowidxs.append(lowidxs[-1]+1+re.search('[atcg]', refam[lowidxs[-1]+1:]).start())
            refas = [refam[lowidxs[1]:lowidxs[2]+1].upper()]
            queryas = [queryam[lowidxs[1]:lowidxs[2]+1].upper()]
            if jbs[jj]<len(chrs1)-1:
                refas.append(refam[lowidxs[3]:lowidxs[4]+1].upper())
                queryas.append(queryam[lowidxs[3]:lowidxs[4]+1].upper())
            sticks = []
            for i in range(lowtimes//2):
                sst = re.search('[ATCG]', queryas[i])
                if sst==None:
                    refas[i], queryas[i] = "", ""
                else:
                    segst = sst.start()
                    segrst = re.search('[ATCG]', queryas[i][::-1]).start()
                    if segrst>0:
                        refas[i], queryas[i] = refas[i][segst:-segrst], queryas[i][segst:-segrst]
                    else:
                        refas[i], queryas[i] = refas[i][segst:], queryas[i][segst:]
                sticks.append("".join(['|' if rc in 'ATCG' and qc in 'ATCG' else ' ' for rc, qc in zip(refas[i], queryas[i])]))

            if jbs[jj]==len(chrs1)-1: # the case of single cut
                si = [0]
            else: # the case of double cut
                si = [ii for ii in range(2) if ends[ii]-starts[ii]>=minseg]
                if len(si)==0: # both mapping segments are short
                    si = [0 if ends[0]-starts[0]>=ends[1]-starts[1] else 1] # use the longer one

            if len(si)==1:
                refaline, queryaline, stickline = refas[si[0]], queryas[si[0]], sticks[si[0]]
            else:
                refaline, queryaline, stickline = f"{refas[0]} {refas[1]}", f"{queryas[0]} {queryas[1]}", f"{sticks[0]} {sticks[1]}"

            _ = fw.write(f">seq{id}\t0\n{queryaline}\n{stickline}\n{refaline}\n")
            _ = ext.write(f">seq{id}\t0\n")
            maps = [re.findall('[ATCG]', querya) for querya in queryas]
            query_pos = 0
            for start, end, mid, map, chrs, cuts, strands, offset, i in zip(starts, ends, mids, maps, [chrs1,chrs2], [cuts1,cuts2], [strands1,strands2], [0,2*refext], [0,1]):
                query_pos += len(mid)
                if i in si:
                    ext.write(f"{genome_file}:>{chrs[jbs[jj]]}")
                    if strands[jbs[jj]]=='-':
                        ext.write(f"_RC\t{ref2len[chrs[jbs[jj]]]-cuts[jbs[jj]]-refext+start-offset}\t")
                    else:
                        ext.write(f"\t{cuts[jbs[jj]]-refext+start-offset}\t")
                    ext.write(f"{query_pos}\t0\n")
                query_pos += len(map)
                if i in si:
                    ext.write(f"{genome_file}:>{chrs[jbs[jj]]}")
                    if strands[jbs[jj]]=='-':
                        ext.write(f"_RC\t{ref2len[chrs[jbs[jj]]]-cuts[jbs[jj]]-refext+end-offset}\t")
                    else:
                        ext.write(f"\t{cuts[jbs[jj]]-refext+end-offset}\t")
                    ext.write(f"{query_pos}\t0\n")

### transform .ext and .alg to .sam
def to_sam(fasta_file):
    # save the genome size to ref2len
    ref2len = {}
    for row in bioframe.read_table('/home/ljw/hg19_with_bowtie2_index/hg19.chrom.sizes', schema=['chrom', 'length']).itertuples():
        ref2len[row.chrom] = row.length
    # scan all .ext and .alg files
    files = [os.path.join(os.path.dirname(fasta_file),file) for file in os.listdir(os.path.dirname(fasta_file)) if file.endswith('.ext') and file.startswith(os.path.basename(fasta_file))]
    with open(f"{fasta_file}.sam", "w") as f:
        for file in files:
            # read all ext results
            with open(file, 'r') as ext:
                exts = ext.readlines()
            # scan alg line by line
            with open(f'{file[:-3]}alg', 'r') as alg:
                # use ei to index records (exts) in .ext file
                ei = -1
                for _, query, mid, ref in more_itertools.batched(alg,4):
                    query_segs = query.strip().split()
                    ref_segs = ref.strip().split()
                    ei += 1
                    ext_segs = exts[ei].strip().split()
                    QNAME = ext_segs[0][1:]
                    MAPIDX = int(ext_segs[1])
                    RNAMEs, strands, ref_starts, ref_ends, query_starts, query_ends, scores = [], [], [], [], [], [], []
                    for i in range(len(query_segs)):
                        ei += 1
                        ext_segs = exts[ei].strip().split()
                        ar = ext_segs[0].split(":")
                        if (ar[-1].endswith("_RC")):
                            RNAMEs.append(ar[-1][1:-3])
                            strands.append("-")
                            ref_ends.append(ref2len[RNAMEs[-1]]-int(ext_segs[1])+1)
                        else:
                            RNAMEs.append(ar[-1][1:])
                            strands.append("+")
                            ref_starts.append(int(ext_segs[1])+1)
                        query_starts.append(int(ext_segs[-2]))
                        scores.append(float(ext_segs[-1]))
                        ei += 1
                        ext_segs = exts[ei].strip().split()
                        if (ar[-1].endswith("_RC")):
                            ref_starts.append(ref2len[RNAMEs[-1]]-int(ext_segs[1])+1)
                        else:
                            ref_ends.append(int(ext_segs[1])+1)
                        query_ends.append(int(ext_segs[-2]))
                        scores[-1] = float(ext_segs[-1]) - scores[-1]
                    # set CIGARs
                    start = 0
                    CIGARs = []
                    for i in range(len(query_segs)):
                        ql = len(query_segs[i])
                        mid_seg = mid[start:start+len(query_segs[i])]
                        query_segs[i] = query_segs[i].lstrip("-")
                        mid_seg = mid_seg[len(mid_seg)-len(query_segs[i]):]
                        ref_segs[i] = ref_segs[i][len(ref_segs[i])-len(query_segs[i]):]
                        query_segs[i] = query_segs[i].rstrip("-")
                        if len(query_segs[i])<len(mid_seg):
                            mid_seg = mid_seg[:len(query_segs[i])-len(mid_seg)]
                            ref_segs[i] = ref_segs[i][:len(query_segs[i])-len(ref_segs[i])]
                        typepre="?"
                        cumsum=0
                        CIGARs.append("")
                        for j in range(len(query_segs[i])):
                            if mid_seg[j]=="|":
                                type="="
                            else:
                                if query_segs[i][j]=="-":
                                    type="D"
                                else:
                                    if ref_segs[i][j]=="-":
                                        type="I"
                                    else:
                                        type="X"
                            if type != typepre:
                                if cumsum > 0:
                                    CIGARs[-1] += f"{cumsum}{typepre}" 
                                    cumsum = 0
                                typepre=type
                            cumsum += 1
                        CIGARs[-1] += f"{cumsum}{type}"
                        start += ql+1
                    # set FLAGs
                    FLAGs = []
                    for i in range(len(query_segs)):
                        FLAGs.append(2)
                        if len(query_segs)>1:
                            FLAGs[-1] += 1
                        if strands[i]=="-":
                            FLAGs[-1] += 16
                        if strands[(i+1)%len(query_segs)]=="-":
                            FLAGs[-1] += 32
                        if i!=0:
                            FLAGs[-1] += 2048
                    # set TAGs
                    TAGs = []
                    for i in range(len(query_segs)):
                        TAGs.append(f"XX:i:{i}\tXZ:i:{MAPIDX}")
                    # set TLENs
                    TLENs = []
                    if len(query_segs)==1:
                        TLENs.append(0)
                    else:
                        unify_chrom = True
                        for i in range(1,len(query_segs)):
                            if RNAMEs[i]!=RNAMEs[0]:
                                for j in range(len(query_segs)):
                                    TLENs.append(0)
                                unify_chrom = False
                        if unify_chrom:
                            leftmost = ref_starts[0]
                            rightmost = ref_ends[0]
                            rightidx = 0
                            for i in range(1,len(query_segs)):
                                if leftmost > ref_starts[i]:
                                    leftmost = ref_starts[i]
                                if rightmost < ref_ends[i]:
                                    rightmost = ref_ends[i]
                                    rightidx = i
                            TLEN = rightmost - leftmost
                            for i in range(len(query_segs)):
                                TLENs.append(TLEN)
                            TLENs[rightidx] = -TLEN
                    # set SEQs
                    SEQs = []
                    for i in range(len(query_segs)):
                        SEQs.append(query_segs[i].replace("-","").upper())
                        if strands[i]=="-":
                            SEQs[-1] = Bio.Seq.Seq(SEQs[-1]).reverse_complement().__str__()
                    # output
                    for i in range(len(query_segs)):
                        f.write(f'{QNAME}\t{FLAGs[i]}\t{RNAMEs[i]}\t{ref_starts[i]}\t255\t{CIGARs[i]}\t{RNAMEs[(i+1)%len(query_segs)]}\t{ref_starts[(i+1)%len(query_segs)]}\t{TLENs[i]}\t{SEQs[i]}\t*\t{TAGs[i]}\n')

def filter_bam(bamfile, region, ori, ord='forward', gap=100, segnums=[2]):
    # get reads in region
    chroms, starts, ends, names, scores, strands, XXs = [], [], [], [], [], [], []
    with pysam.AlignmentFile(bamfile) as sf:
        for read in sf.fetch(region[0], region[1], region[2]):
            chroms.append(read.reference_name)
            starts.append(read.reference_start)
            ends.append(read.reference_end)
            names.append(read.query_name)
            scores.append(read.mapping_quality)
            if read.is_forward:
                strands.append("+")
            else:
                strands.append("-")
            XXs.append(read.get_tag('XX'))
    blocks = pandas.DataFrame({'chrom' : chroms, 'start' : starts, 'end' : ends, 'name' : names, 'score' : scores, 'strand' : strands, 'XX' : XXs})
    blocks = blocks.sort_values(by = ['name', 'XX']).reset_index(drop=True)
    # filter group by ori, ord, gap, segnums
    filter_names = []
    for name, group in blocks.groupby(['name']): 
        if len(group) not in segnums:
            filter_names.append(name)
            continue
        filter = False
        for i in range(len(group)):
            if group['strand'].iloc[i] not in ori[i]:
                filter = True
                break
            if i>0:
                # filter group not consecutive
                if group['XX'].iloc[i]-group['XX'].iloc[i-1]!=1:
                    filter = True
                    break
                # filter gap
                if (ord=='forward' and group['start'].iloc[i]-group['end'].iloc[i-1]<=gap) or (ord=='reverse' and group['start'].iloc[i-1]-group['end'].iloc[i]<=gap):
                    filter = True
                    break
        if filter:
            filter_names.append(name)
    return blocks[~blocks['name'].isin(filter_names)].reset_index(drop=True)



def get_key(isfirst, strand, mapstart, mapend, refstart, refend, mapkeys):
    if isfirst:
        if strand=="+":
            mapkeys[0] = mapend
            maprefkey = mapkeys[0] - refstart
        else:
            mapkeys[0] = mapstart
            maprefkey = refend - mapkeys[0]
    else:
        if strand=="+":
            mapkeys[1] = mapstart
            maprefkey = mapkeys[1] - refstart
        else:
            mapkeys[1] = mapend
            maprefkey = refend - mapkeys[1]
    return maprefkey


def get_indel(fasta_file, count_file, genome_file, region1, strand1, cut1, region2, strand2, cut2, thres=100000, align_mode="up"):
    ref2len = {}
    for row in bioframe.read_table('/home/ljw/hg19_with_bowtie2_index/hg19.chrom.sizes', schema=['chrom', 'length']).itertuples():
        ref2len[row.chrom] = row.length
    files = [os.path.join(os.path.dirname(fasta_file),file) for file in os.listdir(os.path.dirname(fasta_file)) if file.endswith('.ext') and file.startswith(os.path.basename(fasta_file))]
    genome = bioframe.load_fasta(genome_file)
    ref1 = genome[region1[0]].ff.fetch(region1[0],region1[1],region1[2]).upper()
    ref2 = genome[region2[0]].ff.fetch(region2[0],region2[1],region2[2]).upper()
    if strand2=='-':
        ref2 = Bio.Seq.Seq(ref2).reverse_complement().__str__()
    with open(fasta_file, 'r') as fa:
        temp = fa.readlines()
    qnames = [line.rstrip().lstrip('>') for line in temp[::2]]
    temp = bioframe.read_table(count_file, schema=['query','count'])
    queriesfull = temp['query'].values
    counts = temp['count'].values
    indel_file = f"{fasta_file}.{region1[0]}.{region2[0]}.{cut1}.{cut2}.{strand1}.{strand2}.indel"
    with open(indel_file, 'w') as f:
        f.write(f"qname\tcount\tmapidx\tref1\tkey1\tstrand1\tref2\tkey2\tstrand2\tchrom1\tcut1\tstrand1req\tchrom2\tcut2\tstrand2req\tdel1\tins1\tdel2\tins2\tmid_seg\tmmej_down\tmmej_up\ttarget_shift\tref1_seg\tquery1_seg\tmapstart1\tmapend1\tref2_seg\tquery2_seg\tmapstart2\tmapend2\n")
        for file in files:
            # read align results
            with open(file, 'r') as ext:
                exts = ext.readlines()
            ei = -1
            # get indel
            with open(f'{file[:-3]}alg', 'r') as alg:
                qstart = 0
                for _, query, _, ref in more_itertools.batched(alg,4):
                    query_segs = query.strip().split()
                    ref_segs = ref.strip().split()
                    ei += 1
                    ext_segs = exts[ei].strip().split()
                    QNAME = ext_segs[0][1:]
                    MAPIDX = int(ext_segs[1])
                    qstart = qnames.index(QNAME, qstart)
                    count = counts[qstart]
                    queryfull = queriesfull[qstart]
                    RNAMEs, strands, mapends, mapstarts, query_starts, query_ends, scores = [], [], [], [], [], [], [] 
                    for i in range(len(query_segs)):
                        ei += 1
                        ext_segs = exts[ei].split()
                        ar = ext_segs[0].split(":")
                        if ar[-1].endswith("_RC"):
                            RNAMEs.append(ar[-1][1:-3])
                            strands.append("-")
                            mapends.append(ref2len[RNAMEs[-1]]-int(ext_segs[1]))
                        else:
                            RNAMEs.append(ar[-1][1:])
                            strands.append("+")
                            mapstarts.append(int(ext_segs[1]))
                        query_starts.append(int(ext_segs[-2]))
                        scores.append(float(ext_segs[-1]))
                        ei += 1
                        ext_segs = exts[ei].split()
                        if ar[-1].endswith("_RC"):
                            mapstarts.append(ref2len[RNAMEs[-1]]-int(ext_segs[1]))
                        else:
                            mapends.append(int(ext_segs[1]))
                        query_ends.append(int(ext_segs[-2]))
                        scores[-1]=float(ext_segs[-1])-scores[-1]

                    mapkeys = [None,None]
                    map1refkey = get_key(True, strands[0], mapstarts[0], mapends[0], region1[1], region1[2], mapkeys)
                    if len(query_segs)>1:
                        map2refkey = get_key(False, strands[-1], mapstarts[-1], mapends[-1], region2[1], region2[2], mapkeys)
                    else:
                        mapkeys[1]="*"

                    cut1ref = cut1 - region1[1]
                    cut2ref = cut2 - region2[1]
                    if len(query_segs)<2 or RNAMEs[0]!=region1[0] or RNAMEs[-1]!=region2[0] or strands[0]!=strand1 or strands[-1]!=strand2 or abs(map1refkey-cut1ref)>thres or abs(map2refkey-cut2ref)>thres:
                        indels = ["*"] * 5
                        downup = ["*", "*"]
                        target_shift="*"
                    else:
                        if query_starts[-1] > query_ends[0]:
                            downup = [0,0]
                        else:
                            i = 0
                            while i < min(map1refkey,map2refkey):
                                if ref1[map1refkey-i-1]!=ref2[map2refkey-i-1]:
                                    break
                                i += 1
                            downup = [-i]
                            i = 0
                            while i < min(len(ref1)-map1refkey,len(ref2)-map2refkey):
                                if ref1[map1refkey+i]!=ref2[map2refkey+i]:
                                    break
                                i += 1
                            downup.append(i)
                        # get indels
                        indels = [None]*5
                        indels[4] = queryfull[query_ends[0]:query_starts[-1]]
                        if align_mode=="up":
                            target_shift = cut1ref - map1refkey
                        elif align_mode=="down":
                            target_shift = cut2ref - map2refkey
                        else:
                            target_shift = (cut1ref+cut2ref-map1refkey-map2refkey)//2
                        target_shift = max(target_shift, downup[0])
                        target_shift = min(target_shift, downup[1])
                        map1refkey += target_shift
                        map2refkey += target_shift
                        indels[0] = max(cut1ref-map1refkey,0)
                        indels[1] = max(map1refkey-cut1ref,0)
                        indels[2] = max(map2refkey-cut2ref,0)
                        indels[3] = max(cut2ref-map2refkey,0)

                    if len(query_segs)<2:
                        strands.append("*")
                        RNAMEs.append("*")
                        ref_segs.append("*")
                        query_segs.append("*")
                        mapstarts.append("*")
                        mapends.append("*")
                    f.write(f"{QNAME}\t{count}\t{MAPIDX}\t{RNAMEs[0]}\t{mapkeys[0]}\t{strands[0]}\t{RNAMEs[-1]}\t{mapkeys[-1]}\t{strands[-1]}\t{region1[0]}\t{cut1}\t{strand1}\t{region2[0]}\t{cut2}\t{strand2}\t{indels[0]}\t{indels[1]}\t{indels[2]}\t{indels[3]}\t{indels[4]}\t{downup[0]}\t{downup[1]}\t{target_shift}\t{ref_segs[0]}\t{query_segs[0]}\t{mapstarts[0]}\t{mapends[0]}\t{ref_segs[-1]}\t{query_segs[-1]}\t{mapstarts[-1]}\t{mapends[-1]}\n")

    return indel_file

def best_indel(indel_files, out_file):
    df = pandas.concat([pandas.read_csv(indel_file, sep='\t') for indel_file in indel_files]).reset_index(drop=True)
    with open(out_file, 'w') as f:
        f.write("\t".join(df.columns)+"\n")
        for _, group in df.groupby('qname'):
            best_values = [2,0,0]
            for line in group.itertuples():
                ins2 = numpy.inf if line.ins2=="*" else int(line.ins2)
                del2 = numpy.inf if line.del2=="*" else int(line.del2)
                values = [int(line.del1=="*"),ins2,del2]
                if values[0]<best_values[0] or values[0]==best_values[0] and values[1]<best_values[1] or values[0]==best_values[0] and values[1]==best_values[1] and values[2]<best_values[2]:
                    best_values = values
                    best_line = line
            f.write(f"{best_line[1]}")
            for i in range(2,len(best_line)):
                f.write(f"\t{best_line[i]}")
            f.write("\n")
        
def longest_common_substring_cumulative(outfile, n=numpy.arange(9,101), M=2, p=0.25):
    n = n.astype("float64")
    b = M - 2*numpy.log(n)/numpy.log(1/p)
    proupper = (1+p)/(1-p) * (1 - 2/n*numpy.log(n)/numpy.log(1/p))**(-2) * p**(-b-1)
    proupper = numpy.minimum(proupper, 1)
    prolower = 1 - p**(b + 1)
    prolower = numpy.maximum(prolower, 0)
    df = pandas.DataFrame({'length' : n, 'lower' : prolower, 'upper' : prolower})
    df.to_csv(path_or_buf=outfile, sep='\t', header=True, index=False)
    f, ax = matplotlib.pyplot.subplots()
    ax.plot(n, 1-proupper, c="k", ls='-', marker='.')
    ax.legend('lower')
    ax.set(xlim=[min(n),max(n)], xticks=[min(n),max(n)], xticklabels=[min(n),max(n)], ylim=[0,ax.get_ylim()[1]], yticks=[0,ax.get_ylim()[1]], yticklabels=[0,ax.get_ylim()[1]], xlabel='length', ylabel='lower bound of probability of mmej')
    f.tight_layout()
    f.savefig(f"{outfile}.pdf")
    matplotlib.pyplot.close(f)

def get_target_shift(best_indel_slice, align_mode):
    target_shift_left = best_indel_slice['cut1'] - best_indel_slice['key1']
    target_shift_left[best_indel_slice['strand1']=='-'] *= -1
    target_shift_right = best_indel_slice['cut2'] - best_indel_slice['key2'].astype("int64")
    target_shift_left[best_indel_slice['strand2']=='-'] *= -1
    if align_mode=='align_left':
        target_shift = target_shift_left
    elif align_mode=='align_right':
        target_shift = target_shift_right
    else:
        target_shift = (target_shift_left + target_shift_right) // 2
    target_shift = numpy.maximum(target_shift, best_indel_slice['mmej_down'].astype("int64"))
    target_shift = numpy.minimum(target_shift, best_indel_slice['mmej_up'].astype("int64"))
    map1ranges = best_indel_slice[['mapstart1','mapend1']]
    map1ranges.loc[best_indel_slice['strand1']=='+','mapend1'] += target_shift[best_indel_slice['strand1']=='+']
    map1ranges.loc[best_indel_slice['strand1']=='-','mapstart1'] -= target_shift[best_indel_slice['strand1']=='-']
    map2ranges = best_indel_slice[['mapstart2','mapend2']].astype("int64")
    map2ranges.loc[best_indel_slice['strand2']=='+','mapstart2'] += target_shift[best_indel_slice['strand2']=='+']
    map2ranges.loc[best_indel_slice['strand2']=='-','mapend2'] -= target_shift[best_indel_slice['strand2']=='-']
    return target_shift, map1ranges.values, map2ranges.values

def fill_coverage(mapranges, counts, strand):
    mapranges = numpy.array(mapranges)
    counts = numpy.array(counts)
    edges = numpy.concatenate([mapranges[:,0],mapranges[:,1]])
    uni_edges, IC = numpy.unique(edges,return_inverse=True)
    delta_heights = numpy.bincount(IC, numpy.concatenate([counts,-counts]))
    y = numpy.roll(numpy.repeat(numpy.cumsum(delta_heights),2),1)
    x = numpy.repeat(uni_edges,2)
    if strand=='-':
        y = -y
    return x, y

def draw_coverage(best_file, chrs1, cuts1, strands1, chrs2, cuts2, strands2, dis_range, align_mode='align_left'):
    cut12 = list(set(cuts1 + cuts2))
    best_indel = pandas.read_csv(best_file, sep='\t')
    best_indel = best_indel[best_indel['del1']!="*"].reset_index(drop=True)
    mapranges_plus, mapranges_minus, counts_plus, counts_minus = [], [], [], []
    for cc in range(len(chrs1)):
        bool = (best_indel['ref1']==chrs1[cc]) & (best_indel['cut1']==cuts1[cc]) & (best_indel['strand1']==strands1[cc]) & (best_indel['ref2']==chrs2[cc]) & (best_indel['cut2']==cuts2[cc]) & (best_indel['strand2']==strands2[cc])
        best_indel_slice = best_indel[bool].reset_index(drop=True)
        _ , map1ranges , map2ranges = get_target_shift(best_indel_slice, align_mode)
        for strand, mapranges in zip([strands1[cc],strands2[cc]], [map1ranges ,map2ranges]):
            if strand=='+':
                mapranges_plus.extend(mapranges)
                counts_plus.extend(best_indel_slice['count'])
            else:
                mapranges_minus.extend(mapranges)
                counts_minus.extend(best_indel_slice['count'])
    f, ax = matplotlib.pyplot.subplots()
    if mapranges_plus:
        [fx,fy] = fill_coverage(mapranges_plus, counts_plus, '+')
        ax.fill(fx, fy, fc='k', ec=None)
    if mapranges_minus:
        [fx,fy] = fill_coverage(mapranges_minus, counts_minus, '-')
        ax.fill(fx, fy, fc='k', ec=None)
    ax.plot(dis_range, [0,0], c='k', ls='-')
    ylim = ax.get_ylim()
    for cut in cut12:
        ax.plot([cut,cut], ylim, c='k', ls='--')
    ax.set(xlim=dis_range, xticks=[dis_range[0]]+cut12+[dis_range[1]], xticklabels=[dis_range[0]]+cut12+[dis_range[1]], ylim=ylim, yticks=ylim, yticklabels=ylim, title='coverage')
    f.tight_layout()
    f.savefig(f"{best_file}.coverage.pdf")

def get_bool_rv(cut1, strand1, map1ranges, cut2, strand2, map2ranges, cut12):
    if cut1<cut2 or cut1==cut2 and strand1=='+':
        bool_rv = map2ranges[:,0]<map1ranges[:,1]
    else:
        bool_rv = map1ranges[:,0]<map2ranges[:,1]
    if strand1=='+' and (cut1==cut2 and cut1==cut12[0] or cut1>cut2 and strand2=='+'):
        bool_rv = bool_rv | (map2ranges[:,1]>cut12[1])
    if strand1=='-' and (cut1==cut2 and cut1==cut12[1] or cut1<cut2 and strand2=='-'):
        bool_rv = bool_rv | (map2ranges[:,0]<cut12[0])
    
    return bool_rv

def get_title(cut1, strand1, cut2, strand2):
    if strand1=='+' and strand2=='-':
        return 'upstream_inversion'
    else:
        if strand1=='-' and strand2=='+':
            return 'downstream_inversion'
        else:
            if cut1==cut2:
                return 'single_cut'
            else:
                if strand1=='+' and strand2=='+' and cut1 < cut2 or strand1=='-' and strand2=='-' and cut1>cut2:
                    return 'deletion'
                else:
                    return 'duplication'

def sort_by_second_resection(best_indel, align_mode, strand2):
    target_shift, _, _ = get_target_shift(best_indel, align_mode)
    second_resection = best_indel['del2'].astype("int64").values - best_indel['ins2'].astype("int64").values - best_indel['target_shift'].astype("int64").values + target_shift.values
    if strand2=='-':
        idx = numpy.argsort(second_resection)
    else:
        idx = numpy.argsort(-second_resection)
    second_resection = second_resection[idx]
    best_indel = best_indel.iloc[idx].reset_index(drop=True)

    return second_resection, best_indel

def fill_map(mapranges, counts, strand):
    counts = numpy.array(counts)
    y = numpy.cumsum(counts)
    y = numpy.concatenate([[0],numpy.repeat(y[:-1],2),[y[-1]]])
    y = numpy.concatenate([y,numpy.flip(y)])
    x = numpy.concatenate([numpy.repeat(mapranges[:,0],2),numpy.repeat(numpy.flip(mapranges[:,1]),2)])
    if strand=='-':
        y = -y
    return x, y

def draw_large_right_resection(best_file, chr1, cut1, strand1, chr2, cut2, strand2, dis_range, cut12, second_resection_range=[0,99], align_mode='align_left'):
    best_indel = pandas.read_csv(best_file, sep='\t')
    best_indel = best_indel[best_indel['del1']!="*"].reset_index(drop=True)
    bool = (best_indel['ref1']==chr1) & (best_indel['cut1']==cut1) & (best_indel['strand1']==strand1) & (best_indel['ref2']==chr2) & (best_indel['cut2']==cut2) & (best_indel['strand2']==strand2) & (best_indel['mapend2'].astype("int64")>dis_range[0]) & (best_indel['mapstart2'].astype("int64")<dis_range[1])
    best_indel = best_indel[bool].reset_index(drop=True)
    second_resection, best_indel = sort_by_second_resection(best_indel, align_mode, strand2)
    _, map1ranges, map2ranges = get_target_shift(best_indel, align_mode)

    bool_rv = get_bool_rv(cut1, strand1, map1ranges, cut2, strand2, map2ranges, cut12)
    if len(bool_rv)==0:
        warnings.warn("empty reads for this indel")
        return
    bool_rv = bool_rv | (second_resection_range[0]>second_resection) | (second_resection_range[1]<second_resection)
    if not any(~bool_rv):
        return
    counts = best_indel.loc[~bool_rv,'count']

    f, ax = matplotlib.pyplot.subplots()
    for mapranges, strand in zip([map1ranges,map2ranges],[strand1,strand2]):  
        [fx, fy] = fill_map(mapranges[~bool_rv,:], counts, strand)
        ax.fill(fx, fy, fc='k', ec=None)
    disp_total = sum(counts)
    if strand1=='+' and strand2=='+':
        ylim = [0,disp_total]
    elif strand1=='-' and strand2=='-':
        ylim = [-disp_total,0]
    else:
        ylim = [-disp_total,disp_total]
        ax.plot(dis_range, [0,0], c='k', ls='-')
    ax.set(xlim=dis_range, xticks=[dis_range[0],cut12[0],cut12[1],dis_range[1]], xticklabels=[dis_range[0],cut12[0],cut12[1],dis_range[1]], ylim=ylim, yticks=list(set(ylim+[0])), yticklabels=list(set(ylim+[0])), title=get_title(cut1, strand1, cut2, strand2))
    for cut in cut12:
        ax.plot([cut,cut], ylim, c='k', ls='--')
    f.tight_layout()
    f.savefig(f"{best_file}.{chr1}.{chr2}.{cut1}.{cut2}.{strand1}.{strand2}.{second_resection_range[0]}.{second_resection_range[1]}.resection.pdf")


def draw_segs(best_file, chr1, cut1, strand1, chr2, cut2, strand2, dis_range, cut12, align_mode='align_left'):
    best_indel = pandas.read_csv(best_file, sep='\t')
    best_indel = best_indel[best_indel['del1']!="*"].reset_index(drop=True)
    bool = (best_indel['ref1']==chr1) & (best_indel['cut1']==cut1) & (best_indel['strand1']==strand1) & (best_indel['ref2']==chr2) & (best_indel['cut2']==cut2) & (best_indel['strand2']==strand2) & (best_indel['mapend2'].astype("int64")>dis_range[0]) & (best_indel['mapstart2'].astype("int64")<dis_range[1])
    best_indel = best_indel[bool].reset_index(drop=True)
    _, best_indel = sort_by_second_resection(best_indel, align_mode, strand2)
    _, map1ranges, map2ranges = get_target_shift(best_indel, align_mode)

    bool_rv = get_bool_rv(cut1, strand1, map1ranges, cut2, strand2, map2ranges, cut12)
    if not any(~bool_rv):
        return
    counts = best_indel.loc[~bool_rv,'count']

    f, axs = matplotlib.pyplot.subplots(figsize=(40,40), ncols=1, nrows=3, sharex=True)
    for i, mapranges, strand in zip([0,1],[map1ranges,map2ranges],[strand1,strand2]):
        mapranges = mapranges[~bool_rv,:]
        if i==0 and strand=='+' or i==1 and strand=='-':
            mapranges[:,0] = dis_range[0]
        else:
            mapranges[:,1] = dis_range[1]
        [fx, fy] = fill_coverage(mapranges, counts, strand)
        axs[i].fill(fx, fy, fc='k', ec=None)
        ylim = list(axs[i].get_ylim())
        if strand=='+':
            ylim[0] = 0
        else:
            ylim[1] = 0
        for cut in cut12:
            axs[i].plot([cut,cut], ylim, c='k', ls='--')
        axs[i].set(xlim=dis_range, xticks=[dis_range[0],cut12[0],cut12[1],dis_range[1]], xticklabels=[dis_range[0],cut12[0],cut12[1],dis_range[1]], ylim=ylim, yticks=ylim, yticklabels=ylim, title=f'map{i+1}ext')
        axs[2].fill(fx, fy, fc='k', ec=None)

    ylim = list(axs[2].get_ylim())
    if strand1==strand2:
        if strand=='+':
            ylim[0] = 0
        else:
            ylim[1] = 0
    for cut in cut12:
        axs[2].plot([cut,cut], ylim, c='k', ls='--')
    axs[2].set(xlim=dis_range, xticks=[dis_range[0],cut12[0],cut12[1],dis_range[1]], xticklabels=[dis_range[0],cut12[0],cut12[1],dis_range[1]], ylim=ylim, yticks=list(set(ylim+[0])), yticklabels=list(set(ylim+[0])), title='merge')

    f.tight_layout()
    f.savefig(f"{best_file}.{chr1}.{chr2}.{cut1}.{cut2}.{strand1}.{strand2}.line.pdf")



def count_types(best_file, count_file, chr1, cut1, strand1, chr2, cut2, strand2):
    total_reads = sum(bioframe.read_table(count_file, schema=["seq", "count"])["count"])
    best_indel = pandas.read_csv(best_file, sep='\t')
    if chr1==chr2 and cut1==cut2 and strand1==strand2:
        single_seq_num_special = sum(best_indel.loc[(best_indel["ref2"]=="*") & (best_indel["ref1"]==chr1) & (best_indel["strand1"]==strand1) & (cut1>best_indel["mapstart1"]) & (cut1<best_indel["mapend1"]), "count"])
    else:
        single_seq_num_special = 0
    best_indel = best_indel[best_indel["del1"]!="*"].reset_index(drop=True)
    bool = (best_indel["chrom1"]==chr1) & (best_indel["cut1"]==cut1) & (best_indel["strand1req"]==strand1) & (best_indel["chrom2"]==chr2) & (best_indel["cut2"]==cut2) & (best_indel["strand2req"]==strand2)
    best_indel = best_indel[bool].reset_index(drop=True)

    Del= (best_indel["del1"].astype("int64") + best_indel["del2"].astype("int64")).values
    Ins= (best_indel["ins1"].astype("int64") + best_indel["ins2"].astype("int64")).values
    Mid = numpy.array([len(mid_seg) if isinstance(mid_seg,str) else 0 for mid_seg in best_indel["mid_seg"]])
    indelTable = [sum(best_indel.loc[(Mid+Ins==0) & (Del==0),"count"]),sum(best_indel.loc[(Mid+Ins>0) & (Del==0),"count"]),sum(best_indel.loc[(Mid+Ins==0) & (Del>0),"count"]),sum(best_indel.loc[(Mid+Ins>0) & (Del>0),"count"])]
    indelTable[0]+=single_seq_num_special

    f, axs = matplotlib.pyplot.subplots(figsize=(20,20), nrows=3, ncols=1)
    seaborn.histplot(x=Del, weights=best_indel['count'].values, binrange=[0,100], discrete=True, ax=axs[0])
    axs[0].set(ylabel='deletion counts', title=f"{total_reads}")
    seaborn.histplot(x=Ins, weights=best_indel['count'].values, binrange=[0,100], discrete=True, ax=axs[1])
    axs[1].set(ylabel='insertion counts')
    seaborn.barplot(data=pandas.DataFrame({'WT' : [indelTable[0]], 'INS' : [indelTable[1]], 'DEL' : [indelTable[2]], 'INDEL' : [indelTable[3]]}), ax=axs[2])
    axs[2].set(ylabel='indel counts')
    for p in axs[2].patches: 
        axs[2].text(p.get_x() + p.get_width()/2., p.get_height(), f'{int(p.get_height())}', ha="center", va="bottom")
    f.tight_layout()
    f.savefig(f"{best_file}.{chr1}.{chr2}.{cut1}.{cut2}.{strand1}.{strand2}.count.pdf")

def summary_count(best_files, count_files, chrss1, chrss2, cutss1, cutss2, strandss1, strandss2):
    tables = [[] for _ in range(len(chrss1[0]))]
    for best_file, count_file, chrs1, chrs2, cuts1, cuts2, strands1, strands2 in zip(best_files, count_files, chrss1, chrss2, cutss1, cutss2, strandss1, strandss2):
        total_reads = total_reads = sum(bioframe.read_table(count_file, schema=["seq", "count"])["count"])
        best_indel = pandas.read_csv(best_file, sep='\t')
        bool_filter = (best_indel['ref1']==chrs1[0]) & (best_indel['strand1']==strands1[0]) & (best_indel['ref2']=='*') & (best_indel['mapstart1']<cuts1[0]) & (cuts1[0]<best_indel['mapend1']) | (best_indel['del1']!='*')
        filter_total = sum(best_indel.loc[bool_filter,'count'])
        best_indel = best_indel[best_indel["del1"]!="*"].reset_index(drop=True)
        for ii in range(len(chrs1)):
            tables[ii].append([])
            tables[ii][-1].append(os.path.basename(best_file))
            bool = (best_indel['chrom1']==chrs1[ii]) & (best_indel['cut1']==cuts1[ii]) & (best_indel['strand1']==strands1[ii]) &  (best_indel['chrom2']==chrs2[ii]) & (best_indel['cut2']==cuts2[ii]) & (best_indel['strand2']==strands2[ii])
            Del = best_indel['del1'].astype("int64")[bool].values + best_indel['del2'].astype("int64")[bool].values
            Ins = best_indel['ins1'].astype("int64")[bool].values + best_indel['ins2'].astype("int64")[bool].values
            histDel = numpy.bincount(Del, best_indel['count'][bool].values, minlength=101)
            histIns = numpy.bincount(Ins, best_indel['count'][bool].values, minlength=4)
            tables[ii][-1].extend(histIns[1:4])
            tables[ii][-1].extend(histDel[1:101])
            tables[ii][-1].append(sum(histDel[101:]))
            tables[ii][-1].extend(histIns[1:4]/filter_total)
            tables[ii][-1].extend(histDel[1:101]/filter_total)
            tables[ii][-1].append(sum(histDel[101:])/filter_total)
            tables[ii][-1].append(total_reads)
            tables[ii][-1].append(filter_total)
    for ii in range(len(chrss1[0])):
        with open(os.path.join(os.path.dirname(best_files[0]), f"{get_title(cutss1[0][ii],strandss1[0][ii],cutss2[0][ii],strandss2[0][ii])}.summary"), "w") as f:
            VNs = '\t'.join(['sample', 'insert1bp_count', 'insert2bp_count', 'insert3bp_count'] + [f"delete{i}bp_count" for i in range(1,101)] + ['delete_over_100bp_count', 'insert1bp_percent', 'insert2bp_percent', 'insert3bp_percent'] + [f"delete{i}bp_percent" for i in range(1,101)] + ['delete_over_100bp_percent', 'total_num', 'filter_total'])
            f.write(f"{VNs}\n")
            for row in tables[ii]:
                for item in row:
                    if item!=row[0]:
                        f.write('\t')
                    f.write(f"{item}")
                f.write('\n')


