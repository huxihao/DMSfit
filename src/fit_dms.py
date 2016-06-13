from tools import *
from multiprocessing import Pool

def ss_to_lr(ss):
    "Secondary Structure to Loop Region"
    lr = list(ss)
    last = -1
    for i in range(len(lr)):
        if lr[i] == '(':
            last = i
        elif lr[i] == ')':
            if last != -1 and i-last-1 >= 5:
                for j in range(last+1, i):
                    lr[j] = 'o'
                last = -1
    return ''.join(lr)

def extract_data(outfilename, datafolder, col='1', Format='DMS-example', AConly=True, fileNum=-1, spanLen=10):
    filelist = []
    for afile in os.listdir(datafolder):
        if not afile.endswith('.data'):
            continue
        fname = '.'.join(afile.split('.')[:-1])
        if fname not in filelist:
            filelist.append(fname)
    print len(filelist), 'datasets in the path', datafolder
    ## extract pars data of selected genes
    outfile = open(outfilename, 'w')
    data = []
    for fname in filelist:
        seq_file = datafolder+'/'+fname+'.fa'
        if not os.path.exists(seq_file):
            seq_file = datafolder+'/'+fname+'.dot'
        with open(seq_file, 'r') as tempfile:
            gene = tempfile.readline().strip()
            seq = tempfile.readline().strip()
        gene = gene[1:].strip()
#        seq = seq.replace('U','T') ## Use DNA
        seq = seq.replace('T','U') ## Use RNA
        counts = {}
        with open(datafolder+'/'+fname+'.data', 'r') as tempfile:
            for line in tempfile:
                ele = line.split()
                try:
                    if AConly and seq[int(ele[0])-1] not in ['A','C']:
                        continue
                    if Format == 'DMS-example':
                        from math import log, tanh
                        vals = []
                        for i in col:
                            if i == '1' and i != col[0]: ## Denatured
                                vals.append(log(float(ele[int(i)])+1))
                            elif i == '6': ## SASA
                                vals.append(tanh(float(ele[int(i)])-1))
                            else:
                                vals.append(float(ele[int(i)]))
                    else:
                        vals = [float(ele[int(i)]) for i in col]
                    counts[int(ele[0])-1] = vals
                except:
                    pass
        rpkm = mean_std([counts[i][0] for i in counts])[0]
        data.append((rpkm, gene, seq, counts))
    outfile.write("index,tag,seq,gene")
    for i in xrange(len(col)):
        if i == 0:
            outfile.write(",count")
        else:
            outfile.write(",feature"+str(i))
    outfile.write("\n")
    data.sort(reverse=True)
    cc = 0
    for rpkm, gene, seq, counts in data:
        if cc == fileNum:
            break
        print cc, gene, rpkm, len(seq)
        cc += 1
        ## mseq format
        for i in range(len(seq)): ## move for each position
            tag = "1" ## transcript
            if i < spanLen or i >= len(seq)-spanLen:
                tag = "-" + tag ## marked as negative
            seq_code = ""
            seq_code += seq[i].upper()
            out_gene = ""
            if i == 0:
                out_gene = gene
            if i in counts:
                value = ','.join([str(j) for j in counts[i]])
            else:
                tag = "-2" ## mark missing value
                value = ",".join(["0"]*len(col))
            outfile.write(",".join([str(cc), tag, seq_code, out_gene, value])+"\n")
    outfile.close()

def read_table(datafilename, gcol=4, vcol=5):
    import csv
    data_table = csv.reader(open(datafilename, "rb"), delimiter=',')
    data = {}
    head = data_table.next()
    head[gcol] = "[" + head[gcol] + "]"
    head[vcol] = "<" + head[vcol] + ">"
    print "Read:", "-".join(head),
    index = 0
    gene = ""
    values = []
    for row in data_table:
        if row[gcol] != "":
            index += 1
            if gene != "":
                data[gene] = values
            gene = row[gcol]
            values = []
        if int(row[0]) != index:
            warning("Wrong indexes! %s != %s\n"%(row[0], index))
            exit()
        values.append(row[vcol])
    data[gene] = values
    print "with", len(data), "genes."
    return data

def compare_3m(datapath='../data/DMSseq.invitro/without_structures',
               winSize=4, geneNum=10, column=2):
    extract_data('ProbRNA.csv', datapath, col=column, fileNum=geneNum,  spanLen=int(winSize)/2)
    os.system("Rscript ../src/fit_3m.R")
    show(datapath, True)
    with open('log.txt', 'r') as tempfile:
        for line in tempfile:
            show(line)
    os.remove('log.txt')

def cross_validation(path='../data/DMSseq.invitro', data='ProbRNA', exe='fit_po.R', method='MixPoiLin', winSize=10, geneNum=-1, trainOnly=False):
    data = data+'_cv%s_win%s.csv'%(method, winSize)
    extract_data(data, path+'/without_structures', col=12, fileNum=geneNum, spanLen=winSize/2)
    os.system("Rscript ../src/%s CV %s %s %s"%(exe, method, winSize, data))
    if os.path.exists(data+'.log'):
        with open(data+'.log', 'r') as tempfile:
            for line in tempfile:
                show(line)
        os.remove(data+'.log')
        os.remove(data)

def train_predict(path='../data/DMSseq.invitro', data='ProbRNA', exe='fit_po.R', method='MixPoiLin', winSize=10, geneNum=-1, case='TrainPredict'):
    data = data+'_pred%s_win%s.csv'%(method, winSize)
    if case.startswith('Train'):
        extract_data(data, path+'/without_structures', col=12, fileNum=geneNum, spanLen=winSize/2)
        #extract_data(data, path+'/with_structures', col=12, fileNum=-1, spanLen=winSize/2)
        os.system("Rscript ../src/%s Train %s %s %s"%(exe, method, winSize, data))
    if case.endswith('Predict'):
        extract_data(data, path+'/with_structures', col=12, fileNum=-1, spanLen=winSize/2)
        os.system("Rscript ../src/%s Predict %s %s %s"%(exe, method, winSize, data))
    if os.path.exists(data+'.log'):
        with open(data+'.log', 'r') as tempfile:
            for line in tempfile:
                show(line)
        os.remove(data+'.log')
    return data

def example(column="23", exe='fit_po1.R', method='MixPoiLin', winSize=10):
    show(exe)
    show(method)
    show(column)
    show(winSize)
    data = 'eg_col%s_%s_%s_win%s.csv'%(column, exe, method, winSize)
    if os.path.exists(data):
        print '!!Skip:', data
        return data
    path = '../data/DMS.example'
    extract_data(data, path, col=column, fileNum=-1, spanLen=winSize/2)
    os.system("Rscript ../src/%s Train %s %s %s"%(exe, method, winSize, data))
    os.system("Rscript ../src/%s Predict %s %s %s"%(exe, method, winSize, data))
    if os.path.exists(data+'.log'):
        with open(data+'.log', 'r') as tempfile:
            for line in tempfile:
                show(line)
        os.remove(data+'.log')
    return data

def compare_pdb(method):
    mul = []
    path = '../data/DMS.example'
    from cmp_ss import get_value
    for dname in os.listdir(path):
        if not dname.endswith('.data'):
            continue
        if not dname.startswith('yeast.rRNA_RDN18'):
            continue
#        if not dname.startswith('25S.yeast.domain3'):
#            continue
        one = get_value(path, dname[:-5], method=method[:-4], chs='AC')
        mul += one
    out = list(zip(*mul))
    for pred in out[2:]:
        pair = [(ss,pd) for ss, sa, pd in zip(out[0], out[1], pred) if sa>1 or ss==False]
        r,p = zip(*pair)
        area, x,y,z = performance(r, p)
        show(area)
        #show(correlation(out[1], pred, rank=True))
    show()

def fig1(m='nb'):
    show('Figure1\n')
    for k in [10]:
        example(column=12, exe='fit_'+m+'1.R', method='MixPoiLin', winSize=k)
        example(column=13, exe='fit_'+m+'1.R', method='MixPoiLin', winSize=k)
        example(column=12, exe='fit_'+m+'2.R', method='MixPoiLin', winSize=k)
        example(column=13, exe='fit_'+m+'2.R', method='MixPoiLin', winSize=k)

def fig2():
    show('Figure2\n')
    for m in ['po1', 'po2', 'nb1', 'nb2']:
        for k in [2, 4, 10, 20]:
            compare_pdb(example(column="2", exe='fit_'+m+'.R', method='MPL', winSize=k))
            compare_pdb(example(column="21", exe='fit_'+m+'.R', method='MPL', winSize=k))
            compare_pdb(example(column="3", exe='fit_'+m+'.R', method='MPL', winSize=k))
            compare_pdb(example(column="31", exe='fit_'+m+'.R', method='MPL', winSize=k))
#        compare_pdb(example(column="24", exe='fit_'+m+'.R', method='MPL', winSize=k))
#        compare_pdb(example(column="26", exe='fit_'+m+'.R', method='MPL', winSize=k))
#        compare_pdb(example(column="246", exe='fit_'+m+'.R', method='MPL', winSize=k))

def main(para):
#    fig1()
    fig2()

def others(para):
    for k in [10]:
        example(exe='fit_po1.R', method='MPLStepwise0', winSize=k)
#        example(exe='fit_nb2.R', method='MixPoiLin', winSize=k)
#        cross_validation(path='../data/DMSseq.invitro', data='DMS_vitro', exe='fit_nb1.R', method='MPLStepwise1', winSize=k, geneNum=100)
#        cross_validation(path='../data/DMSseq.invivo', data='DMS_vivo', exe='fit_nb1.R', method='MixPoiLin', winSize=k, geneNum=100)

#        train_predict(path='../data/DMSseq.invitro', data='DMS_vitro_nb1', exe='fit_nb1.R', method='MixPoiLin', winSize=k, geneNum=100, case='TrainPredict')
#        train_predict(path='../data/DMSseq.invitro', data='DMS_vitro_nb1', exe='fit_nb1.R', method='MPLStepwise1', winSize=k, geneNum=100)
#        train_predict(path='../data/DMSseq.invitro', data='DMS_vitro_nb1', exe='fit_nb1.R', method='MPLStepwise2', winSize=k, geneNum=100)
#        train_predict(path='../data/DMSseq.invitro', data='DMS_vitro_nb1', exe='fit_nb1.R', method='MPLStepwise3', winSize=k, geneNum=100)

#        train_predict(path='../data/DMSseq.invivo', data='DMS_vivo_nb1', exe='fit_nb1.R', method='MixPoiLin', winSize=k, geneNum=100, case='TrainPredict')
#        train_predict(path='../data/DMSseq.invivo', data='DMS_vivo_nb1', exe='fit_nb1.R', method='MPLStepwise1', winSize=k, geneNum=100)
#        train_predict(path='../data/DMSseq.invivo', data='DMS_vivo_nb1', exe='fit_nb1.R', method='MPLStepwise2', winSize=k, geneNum=100)
#        train_predict(path='../data/DMSseq.invivo', data='DMS_vivo_nb1', exe='fit_nb1.R', method='MPLStepwise3', winSize=k, geneNum=100)
    show('DMSseq_control')
    compare_3m('../data/DMSseq.invitro/without_structures', winSize=20, geneNum=100, column=1)
    show('DMSseq_vitro')
    compare_3m('../data/DMSseq.invitro/without_structures', winSize=20, geneNum=100, column=2)
    show('DMSseq_vivo')
    compare_3m('../data/DMSseq.invivo/without_structures', winSize=20, geneNum=100, column=2)
    for k in [2, 10, 20]:
        train_predict(path='../data/DMSseq.invitro', data='DMS_control', winSize=k, geneNum=100, column=1)
        train_predict(path='../data/DMSseq.invitro', data='DMS_vitro', winSize=k, geneNum=100, column=2)
        train_predict(path='../data/DMSseq.invivo', data='DMS_vivo', winSize=k, geneNum=100, column=2)

if __name__ == "__main__": main_fun(main)
