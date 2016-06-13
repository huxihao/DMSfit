from tools import *

def get_data1(fil, col):
    infile = open(fil, 'r')
    data = []
    for line in infile:
        ele = line.split('\t')
        data.append(ele[col].strip())
    infile.close()
    return data

def get_data2(fil, nam, col=1):
    infile = open(fil, 'r')
    ix = -1; cc = 0
    line = infile.readline()
    ele = line.strip().split(',')
    for i in xrange(len(ele)):
        if ele[i] == col:
            col = i
            break
    print 'Read', nam, 'in column', col
    name = ''
    out = []
    for line in infile:
        ele = line.split(',')
        idx = int(ele[0])
        if ix == idx:
            cc += 1
        else:
            ix = idx
            cc = 0
            Name = ele[3].strip()
            if Name != '':
                name = Name
        if name == nam:
            out.append(float(ele[col]))
    infile.close()
    return out

def get_value(path, dname, method='nb1_MPL_win60', chs='ACGU', case='rRNA'):
    dfile = path+'/'+dname+'.data'
    idx = get_data1(dfile, col=0)
    ss = get_data1(dfile, col=4)
    nu = get_data1(dfile, col=5)
    sasa = get_data1(dfile, col=6)
    dms_co = get_data1(dfile, col=1)
    dms_vt = get_data1(dfile, col=2)
    dms_vv = get_data1(dfile, col=3)
    real = [i=='0' for i,j in zip(ss, nu) if j in chs] 
    solv = [float(i) for i,j in zip(sasa, nu) if j in chs]
    co = [float(i) for i,j in zip(dms_co, nu) if j in chs]
    vt = [float(i) for i,j in zip(dms_vt, nu) if j in chs]
    vv = [float(i) for i,j in zip(dms_vv, nu) if j in chs]
    if not os.path.exists(method+'.csv'):
        val_co = [float(i) for i in dms_co]
        val_vt = [float(i) for i in dms_vt]
        val_vv = [float(i) for i in dms_vv]
        top95_co = sorted(val_co)[int(0.95*len(nu))]
        top95_vt = sorted(val_vt)[int(0.95*len(nu))]
        top95_vv = sorted(val_vv)[int(0.95*len(nu))]
        pos = range(len(nu))
        if case == 'rRNA':
            co = [min(top95_co, val_co[i]) for i in pos if nu[i] in chs]
            vt = [min(top95_vt, val_vt[i]) for i in pos if nu[i] in chs]
            vv = [min(top95_vv, val_vv[i]) for i in pos if nu[i] in chs]
        elif case == 'mRNA':
            co = [val_co[i]/max([val_co[j] for j in xrange(max(0,i-25), min(i+25, len(nu)))]) for i in pos if nu[i] in chs]
            vt = [val_vt[i]/max([val_vt[j] for j in xrange(max(0,i-25), min(i+25, len(nu)))]) for i in pos if nu[i] in chs]
            vv = [val_vv[i]/max([val_vv[j] for j in xrange(max(0,i-25), min(i+25, len(nu)))]) for i in pos if nu[i] in chs]
        return zip(real, solv, co, vt, vv)
    mpl_co = get_data2(method+'.csv', dname, col="prob")
    sco = [mpl_co[int(i)-1] for i in idx]
    mco = [i for i,j in zip(sco, nu) if j in chs]
    return zip(real, solv, co, vt, vv, mco)

def get_value2(path, dname, method='nb1_predMPLStepwise_win60', chs='ACGU', case='rRNA'):
    dfile = path+'/'+dname+'.data'
    idx = get_data1(dfile, col=0)
    ss = get_data1(dfile, col=4)
    nu = get_data1(dfile, col=5)
    sasa = get_data1(dfile, col=6)
    dms_co = get_data1(dfile, col=1)
    dms_vt = get_data1(dfile, col=2)
    dms_vv = get_data1(dfile, col=3)
    real = [i=='0' for i,j in zip(ss, nu) if j in chs] 
    solv = [float(i) for i,j in zip(sasa, nu) if j in chs]
    co = [float(i) for i,j in zip(dms_co, nu) if j in chs]
    vt = [float(i) for i,j in zip(dms_vt, nu) if j in chs]
    vv = [float(i) for i,j in zip(dms_vv, nu) if j in chs]
    if not os.path.exists(method+'.csv'):
        val_co = [float(i) for i in dms_co]
        val_vt = [float(i) for i in dms_vt]
        val_vv = [float(i) for i in dms_vv]
        top95_co = sorted(val_co)[int(0.95*len(nu))]
        top95_vt = sorted(val_vt)[int(0.95*len(nu))]
        top95_vv = sorted(val_vv)[int(0.95*len(nu))]
        pos = range(len(nu))
        if case == 'rRNA':
            co = [min(top95_co, val_co[i]) for i in pos if nu[i] in chs]
            vt = [min(top95_vt, val_vt[i]) for i in pos if nu[i] in chs]
            vv = [min(top95_vv, val_vv[i]) for i in pos if nu[i] in chs]
        elif case == 'mRNA':
            co = [val_co[i]/max([val_co[j] for j in xrange(max(0,i-25), min(i+25, len(nu)))]) for i in pos if nu[i] in chs]
            vt = [val_vt[i]/max([val_vt[j] for j in xrange(max(0,i-25), min(i+25, len(nu)))]) for i in pos if nu[i] in chs]
            vv = [val_vv[i]/max([val_vv[j] for j in xrange(max(0,i-25), min(i+25, len(nu)))]) for i in pos if nu[i] in chs]
        return zip(real, solv, co, vt, vv)
    mpl_co = get_data2(method+'.csv', dname, col='pbv')
    mpl_vt = get_data2(method+'.csv', dname, col='pbs')
    sco = [mpl_co[int(i)-1] for i in idx]
    svt = [mpl_vt[int(i)-1] for i in idx]
    mco = [i for i,j in zip(sco, nu) if j in chs]
    mvt = [i for i,j in zip(svt, nu) if j in chs]
    return zip(real, solv, co, vt, vv, mco, mvt)


def main(para):
    for m in [#'DMS_vitro_nb1_predMixPoiLin_win20',
              #'DMS_vitro_nb1_predMPLStepwise1_win20',
              #'DMS_vitro_nb1_predMPLStepwise2_win20',
              #'DMS_vitro_nb1_predMPLStepwise3_win20',
              #'DMS_vivo_nb1_predMixPoiLin_win20', 
              #'DMS_vivo_nb1_predMPLStepwise1_win20',
              #'DMS_vivo_nb1_predMPLStepwise2_win20',
              #'DMS_vivo_nb1_predMPLStepwise3_win20',
              'example_win10']:
        show(m, True)
        show()
        for c in ['A','C','G','U','AC','ACGU']:
            show(c)
            path = para['DataPath']+'/DMS.example'
            mul = []
            for dname in os.listdir(path):
                if not dname.endswith('.data'):
                    continue
                if not dname.startswith('yeast.rRNA_RDN18'):
                    continue
                one = get_value(path, dname[:-5], method=m, chs=c)
                mul += one
            out = list(zip(*mul))
            for pred in out[2:]:
                pair = [(ss,pd) for ss, sa, pd in zip(out[0], out[1], pred) if sa>2 or ss==False]
                r,p = zip(*pair)
                area, x,y,z = performance(r, p)
                show(area)
                #show(correlation(out[1], pred, rank=True))
            show()

if __name__ == '__main__': main_fun(main)
