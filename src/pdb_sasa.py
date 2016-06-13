from tools import *

## PyMol Setting
import __main__
__main__.pymol_argv = ['pymol', '-qc']
import pymol

## http://pldserver1.biochem.queensu.ca/~rlc/work/teaching/BCHM823/pymol/alignment/
## http://www.biostars.org/p/42474/
from pymol import cmd
## http://pymolwiki.org/index.php/InterfaceResidues
from pymol import stored
def interfaceResidues(cmpx, cA='c. A', cB='c. B', cutoff=1.0, selName="interface"):
    """
    interfaceResidues -- finds 'interface' residues between two chains in a complex.
 
    PARAMS
        cmpx
            The complex containing cA and cB
 
        cA
            The first chain in which we search for residues at an interface
            with cB
 
        cB
            The second chain in which we search for residues at an interface
            with cA
 
        cutoff
            The difference in area OVER which residues are considered
            interface residues.  Residues whose dASA from the complex to
            a single chain is greater than this cutoff are kept.  Zero
            keeps all residues.
 
        selName
            The name of the selection to return.
 
    RETURNS
        * A selection of interface residues is created and named
            depending on what you passed into selName
        * An array of values is returned where each value is:
            ( modelName, residueNumber, dASA )
 
    NOTES
        If you have two chains that are not from the same PDB that you want
        to complex together, use the create command like:
            create myComplex, pdb1WithChainA or pdb2withChainX
        then pass myComplex to this script like:
            interfaceResidues myComlpex, c. A, c. X
 
        This script calculates the area of the complex as a whole.  Then,
        it separates the two chains that you pass in through the arguments
        cA and cB, alone.  Once it has this, it calculates the difference
        and any residues ABOVE the cutoff are called interface residues.
 
    AUTHOR:
        Jason Vertrees, 2009.        
    """
    # Save user's settings, before setting dot_solvent
    oldDS = cmd.get("dot_solvent")
    cmd.set("dot_solvent", 1)
 
    # set some string names for temporary objects/selections
    tempC, selName1 = "tempComplex", selName+"1"
    chA, chB = "chA", "chB"
 
    # operate on a new object & turn off the original
    cmd.create(tempC, cmpx)
    cmd.disable(cmpx)
 
    # remove cruft and inrrelevant chains
    cmd.remove(tempC + " and not (polymer and (%s or %s))" % (cA, cB))
 
    # get the area of the complete complex
    cmd.get_area(tempC, load_b=1)
    # copy the areas from the loaded b to the q, field.
    cmd.alter(tempC, 'q=b')
 
    # extract the two chains and calc. the new area
    # note: the q fields are copied to the new objects
    # chA and chB
    cmd.extract(chA, tempC + " and (" + cA + ")")
    cmd.extract(chB, tempC + " and (" + cB + ")")
    cmd.get_area(chA, load_b=1)
    cmd.get_area(chB, load_b=1)
 
    # update the chain-only objects w/the difference
    cmd.alter( "%s or %s" % (chA,chB), "b=b-q" )
 
    # The calculations are done.  Now, all we need to
    # do is to determine which residues are over the cutoff
    # and save them.
    stored.r, rVal, seen = [], [], []
    cmd.iterate('%s or %s' % (chA, chB), 'stored.r.append((model,resi,b))')
 
    cmd.enable(cmpx)
    cmd.select(selName1, None)
    for (model,resi,diff) in stored.r:
        key=resi+"-"+model
        if abs(diff)>=float(cutoff):
            if key in seen: continue
            else: seen.append(key)
            rVal.append( (model,resi,diff) )
            # expand the selection here; I chose to iterate over stored.r instead of
            # creating one large selection b/c if there are too many residues PyMOL
            # might crash on a very large selection.  This is pretty much guaranteed
            # not to kill PyMOL; but, it might take a little longer to run.
            cmd.select( selName1, selName1 + " or (%s and i. %s)" % (model,resi))
 
    # this is how you transfer a selection to another object.
    cmd.select(selName, cmpx + " in " + selName1)
    # clean up after ourselves
    cmd.delete(selName1)
    cmd.delete(chA)
    cmd.delete(chB)
    cmd.delete(tempC)
    # show the selection
    cmd.enable(selName)
 
    # reset users settings
    cmd.set("dot_solvent", oldDS)
 
    return rVal

## Ref: http://www.pymolwiki.org/index.php/Get_Area
def get_SASA(pdbfile, pdb1='pdb1'):
    ## set
    cmd.load(pdbfile, pdb1)
    oldDS = cmd.get("dot_solvent")
    cmd.h_add()
    cmd.flag("ignore", "none")
    cmd.set("dot_solvent", 1)
    cmd.set("dot_density", 2)
    cmd.set("solvent_radius", 3)

    ## calculate
    area = cmd.get_area(pdb1, load_b=1)
    stored.r = []
    cmd.iterate(pdb1, 'stored.r.append((model,chain,resi,resn,name,b))')
    areas = {}
    for model, chain, idx, char, name, sasa in stored.r:
        if idx == '653':
            print chain, idx, char, name, sasa
        base = areas.get((chain,idx,char), 0)
        if (char == 'A' and name == 'N1') or (char == 'C' and name == 'N3'):
            base += float(sasa)
        areas[(chain,idx,char)] = base

    ## reset
    cmd.set("dot_solvent", oldDS)
    cmd.delete(pdb1)
    print pdbfile, area, sum(areas.values())
    return areas

Seq_18S = '''
UAUCUGGUUGAUCCUGCCAGUAGUCAUAUGCUUGUCUCAAAGAUUAAGCCAUGCAUGUCUAAGUAUAAGCAAUUUAUACA
GUGAAACUGCGAAUGGCUCAUUAAAUCAGUUAUCGUUUAUUUGAUAGUUCCUUUACUACAUGGUAUAACUGUGGUAAUUC
UAGAGCUAAUACAUGCUUAAAAUCUCGACCCUUUGGAAGAGAUGUAUUUAUUAGAUAAAAAAUCAAUGUCUUCGGACUCU
UUGAUGAUUCAUAAUAACUUUUCGAAUCGCAUGGCCUUGUGCUGGCGAUGGUUCAUUCAAAUUUCUGCCCUAUCAACUUU
CGAUGGUAGGAUAGUGGCCUACCAUGGUUUCAACGGGUAACGGGGAAUAAGGGUUCGAUUCCGGAGAGGGAGCCUGAGAA
ACGGCUACCACAUCCAAGGAAGGCAGCAGGCGCGCAAAUUACCCAAUCCUAAUUCAGGGAGGUAGUGACAAUAAAUAACG
AUACAGGGCCCAUUCGGGUCUUGUAAUUGGAAUGAGUACAAUGUAAAUACCUUAACGAGGAACAAUUGGAGGGCAAGUCU
GGUGCCAGCAGCCGCGGUAAUUCCAGCUCCAAUAGCGUAUAUUAAAGUUGUUGCAGUUAAAAAGCUCGUAGUUGAACUUU
GGGCCCGGUUGGCCGGUCCGAUUUUUUCGUGUACUGGAUUUCCAACGGGGCCUUUCCUUCUGGCUAACCUUGAGUCCUUG
UGGCUCUUGGCGAACCAGGACUUUUACUUUGAAAAAAUUAGAGUGUUCAAAGCAGGCGUAUUGCUCGAAUAUAUUAGCAU
GGAAUAAUAGAAUAGGACGUUUGGUUCUAUUUUGUUGGUUUCUAGGACCAUCGUAAUGAUUAAUAGGGACGGUCGGGGGC
AUCAGUAUUCAAUUGUCAGAGGUGAAAUUCUUGGAUUUAUUGAAGACUAACUACUGCGAAAGCAUUUGCCAAGGACGUUU
UCAUUAAUCAAGAACGAAAGUUAGGGGAUCGAAGAUGAUCAGAUACCGUCGUAGUCUUAACCAUAAACUAUGCCGACUAG
GGAUCGGGUGGUGUUUUUUUAAUGACCCACUCGGCACCUUACGAGAAAUCAAAGUCUUUGGGUUCUGGGGGGAGUAUGGU
CGCAAGGCUGAAACUUAAAGGAAUUGACGGAAGGGCACCACCAGGAGUGGAGCCUGCGGCUUAAUUUGACUCAACACGGG
GAAACUCACCAGGUCCAGACACAAUAAGGAUUGACAGAUUGAGAGCUCUUUCUUGAUUUUGUGGGUGGUGGUGCAUGGCC
GUUCUUAGUUGGUGGAGUGAUUUGUCUGCUUAAUUGCGAUAACGAACGAGACCUUAACCUACUAAAUAGUGGUGCUAGCA
UUUGCUGGUUAUCCACUUCUUAGAGGGACUAUCGGUUUCAAGCCGAUGGAAGUUUGAGGCAAUAACAGGUCUGUGAUGCC
CUUAGACGUUCUGGGCCGCACGCGCGCUACACUGACGGAGCCAGCGAGUCUAACCUUGGCCGAGAGGUCUUGGUAAUCUU
GUGAAACUCCGUCGUGCUGGGGAUAGAGCAUUGUAAUUAUUGCUCUUCAACGAGGAAUUCCUAGUAAGCGCAAGUCAUCA
GCUUGCGUUGAUUACGUCCCUGCCCUUUGUACACACCGCCCGUCGCUAGUACCGAUUGAAUGGCUUAGUGAGGCCUCAGG
AUCUGCUUAGAGAAGGGGGCAACUCCAUCUCAGAGCGGAGAAUUUGGACAAACUUGGUCAUUUAGAGGAACUAAAAGUCG
UAACAAGGUUUCCGUAGGUGAACCUGCGGAAGGAUCAUUA
'''
Seq_25S = '''
GUUUGACCUCAAAUCAGGUAGGAGUACCCGCUGAACUUAAGCAUAUCAAUAAGCGGAGGAAAAGAAACCAACCGGGAUUG
CCUUAGUAACGGCGAGUGAAGCGGCAAAAGCUCAAAUUUGAAAUCUGGUACCUUCGGUGCCCGAGUUGUAAUUUGGAGAG
GGCAACUUUGGGGCCGUUCCUUGUCUAUGUUCCUUGGAACAGGACGUCAUAGAGGGUGAGAAUCCCGUGUGGCGAGGAGU
GCGGUUCUUUGUAAAGUGCCUUCGAAGAGUCGAGUUGUUUGGGAAUGCAGCUCUAAGUGGGUGGUAAAUUCCAUCUAAAG
CUAAAUAUUGGCGAGAGACCGAUAGCGAACAAGUACAGUGAUGGAAAGAUGAAAAGAACUUUGAAAAGAGAGUGAAAAAG
UACGUGAAAUUGUUGAAAGGGAAGGGCAUUUGAUCAGACAUGGUGUUUUGUGCCCUCUGCUCCUUGUGGGUAGGGGAAUC
UCGCAUUUCACUGGGCCAGCAUCAGUUUUGGUGGCAGGAUAAAUCCAUAGGAAUGUAGCUUGCCUCGGUAAGUAUUAUAG
CCUGUGGGAAUACUGCCAGCUGGGACUGAGGACUGCGACGUAAGUCAAGGAUGCUGGCAUAAUGGUUAUAUGCCGCCCGU
CUUGAAACACGGACCAAGGAGUCUAACGUCUAUGCGAGUGUUUGGGUGUAAAACCCAUACGCGUAAUGAAAGUGAACGUA
GGUUGGGGCCUCGCAAGAGGUGCACAAUCGACCGAUCCUGAUGUCUUCGGAUGGAUUUGAGUAAGAGCAUAGCUGUUGGG
ACCCGAAAGAUGGUGAACUAUGCCUGAAUAGGGUGAAGCCAGAGGAAACUCUGGUGGAGGCUCGUAGCGGUUCUGACGUG
CAAAUCGAUCGUCGAAUUUGGGUAUAGGGGCGAAAGACUAAUCGAACCAUCUAGUAGCUGGUUCCUGCCGAAGUUUCCCU
CAGGAUAGCAGAAGCUCGUAUCAGUUUUAUGAGGUAAAGCGAAUGAUUAGAGGUUCCGGGGUCGAAAUGACCUUGACCUA
UUCUCAAACUUUAAAUAUGUAAGAAGUCCUUGUUACUUAAUUGAACGUGGACAUUUGAAUGAAGAGCUUUUAGUGGGCCA
UUUUUGGUAAGCAGAACUGGCGAUGCGGGAUGAACCGAACGUAGAGUUAAGGUGCCGGAAUACACGCUCAUCAGACACCA
CAAAAGGUGUUAGUUCAUCUAGACAGCCGGACGGUGGCCAUGGAAGUCGGAAUCCGCUAAGGAGUGUGUAACAACUCACC
GGCCGAAUGAACUAGCCCUGAAAAUGGAUGGCGCUCAAGCGUGUUACCUAUACUCUACCGUCAGGGUUGAUAUGAUGCCC
UGACGAGUAGGCAGGCGUGGAGGUCAGUGACGAAGCCUAGACCGUAAGGUCGGGUCGAACGGCCUCUAGUGCAGAUCUUG
GUGGUAGUAGCAAAUAUUCAAAUGAGAACUUUGAAGACUGAAGUGGGGAAAGGUUCCACGUCAACAGCAGUUGGACGUGG
GUUAGUCGAUCCUAAGAGAUGGGGAAGCUCCGUUUCAAAGGCCUGAUUUUAUGCAGGCCACCAUCGAAAGGGAAUCCGGU
UAAGAUUCCGGAACCUGGAUAUGGAUUCUUCACGGUAACGUAACUGAAUGUGGAGACGUCGGCGCGAGCCCUGGGAGGAG
UUAUCUUUUCUUCUUAACAGCUUAUCACCCCGGAAUUGGUUUAUCCGGAGAUGGGGUCUUAUGGCUGGAAGAGGCCAGCA
CCUUUGCUGGCUCCGGUGCGCUUGUGACGGCCCGUGAAAAUCCACAGGAAGGAAUAGUUUUCAUGCCAGGUCGUACUGAU
AACCGCAGCAGGUCUCCAAGGUGAACAGCCUCUAGUUGAUAGAAUAAUGUAGAUAAGGGAAGUCGGCAAAAUAGAUCCGU
AACUUCGGGAUAAGGAUUGGCUCUAAGGGUCGGGUAGUGAGGGCCUUGGUCAGACGCAGCGGGCGUGCUUGUGGACUGCU
UGGUGGGGCUUGCUCUGCUAGGCGGACUACUUGCGUGCCUUGUUGUAGACGGCCUUGGUAGGUCUCUUGUAGACCGUCGC
UUGCUACAAUUAACGAUCAACUUAGAACUGGUACGGACAAGGGGAAUCUGACUGUCUAAUUAAAACAUAGCAUUGCGAUG
GUCAGAAAGUGAUGUUGACGCAAUGUGAUUUCUGCCCAGUGCUCUGAAUGUCAAAGUGAAGAAAUUCAACCAAGCGCGGG
UAAACGGCGGGAGUAACUAUGACUCUCUUAAGGUAGCCAAAUGCCUCGUCAUCUAAUUAGUGACGCGCAUGAAUGGAUUA
ACGAGAUUCCCACUGUCCCUAUCUACUAUCUAGCGAAACCACAGCCAAGGGAACGGGCUUGGCAGAAUCAGCGGGGAAAG
AAGACCCUGUUGAGCUUGACUCUAGUUUGACAUUGUGAAGAGACAUAGAGGGUGUAGAAUAAGUGGGAGCUUCGGCGCCA
GUGAAAUACCACUACCUUUAUAGUUUCUUUACUUAUUCAAUGAAGCGGAGCUGGAAUUCAUUUUCCACGUUCUAGCAUUC
AAGGUCCCAUUCGGGGCUGAUCCGGGUUGAAGACAUUGUCAGGUGGGGAGUUUGGCUGGGGCGGCACAUCUGUUAAACGA
UAACGCAGAUGUCCUAAGGGGGGCUCAUGGAGAACAGAAAUCUCCAGUAGAACAAAAGGGUAAAAGCCCCCUUGAUUUUG
AUUUUCAGUGUGAAUACAAACCAUGAAAGUGUGGCCUAUCGAUCCUUUAGUCCCUCGGAAUUUGAGGCUAGAGGUGCCAG
AAAAGUUACCACAGGGAUAACUGGCUUGUGGCAGUCAAGCGUUCAUAGCGACAUUGCUUUUUGAUUCUUCGAUGUCGGCU
CUUCCUAUCAUACCGAAGCAGAAUUCGGUAAGCGUUGGAUUGUUCACCCACUAAUAGGGAACGUGAGCUGGGUUUAGACC
GUCGUGAGACAGGUUAGUUUUACCCUACUGAUGAAUGUUACCGCAAUAGUAAUUGAACUUAGUACGAGAGGAACAGUUCA
UUCGGAUAAUUGGUUUUUGCGGCUGUCUGAUCAGGCAUUGCCGCGAAGCUACCAUCCGCUGGAUUAUGGCUGAACGCCUC
UAAGUCAGAAUCCAUGCUAGAACGCGGUGAUUUCUUUGCUCCACACAAUAUAGAUGGAUACGAAUAAGGCGUCCUUGUGG
CGUCGCUGAACCAUAGCAGGCUAGCAACGGUGCACUUGGCGGAAAGGCCUUGGGUGCUUGCUGGCGAAUUGCAAUGUCAU
UUUGCGUGGGGAUAAAUCAUUUGUAUACGACUUAGAUGUACAACGGGGUAUUGUAAGCAGUAGAGUAGCCUUGUUGUUAC
GAUCUGCUGAGAUUAAGCCUUUGUUGUCUGAUUUGU
'''

def format_data(data, chain, refseq):
    vals = [0]*len(refseq)
    for ch, idx, nu in data:
        if ch != chain:
            continue
        if nu not in ['A','C','G','U']:
            continue
        i = int(idx)-1
        if nu != refseq[i]:
            print 'Throw', ch, idx, nu, refseq[i]
        vals[i] = data[(ch, idx, nu)]
    return refseq, vals

def update_dms(sasa, path, prefix, suffix='.data'):
    path1 = '/DMSseq.invitro/with_structures'
    path2 = '/DMSSeq.invivo/with_structures'
    path3 = '/DMS.example'
    files = os.listdir(path+path1)
    for f in files:
        if f.startswith(prefix) and f.endswith(suffix):
            data1 = []
            with open(path+path1+'/'+f, 'r') as tempfile:
                for line in tempfile:
                    data1.append(line.strip().split('\t'))
            data2 = []
            with open(path+path2+'/'+f, 'r') as tempfile:
                for line in tempfile:
                    data2.append(line.strip().split('\t'))
            idx, rd1, rd2, ss, nu = zip(*data1)
            id_, rd_, rd3, s_, n_ = zip(*data2)
            if len(idx) != len(id_):
                ## fill 0 in rd3
                i = 0; j = 0
                new_rd3 = []
                while i < len(idx) and j < len(id_):
                    if idx[i] == id_[j]:
                        new_rd3.append(rd3[j])
                        j = j + 1
                    else:
                        print 'Fix', f, idx[i], id_[j]
                        new_rd3.append('0')
                    i = i+1
                rd3 = new_rd3
            Seq, Val = sasa
            seq = ''.join(nu)
            st = Seq.find(seq)
            show([f, len(Seq), len(seq), st], True)
            assert Seq[st:(st+len(seq))] == seq
            sa = [str(i) for i in Val[st:(st+len(seq))]]
            sa1 = [str(int(i>2)) for i in Val[st:(st+len(seq))]]
            pos = [str(i/float(len(Seq))) for i in xrange(st, st+len(seq))]
            with open(path+path3+'/'+f, 'w') as tempfile:
                for e in zip(idx, rd1, rd2, rd3, ss, nu, sa):
                    tempfile.write('\t'.join(e)+'\n')

def main(para):
    pdb1 = get_SASA(para['DataPath']+'/DMS.example/3U5B.pdb')
    sasa_18S = format_data(pdb1, '2', Seq_18S.replace('\n',''))
    update_dms(sasa_18S, para['DataPath'], 'yeast.rRNA_RDN18')

    pdb1 = get_SASA(para['DataPath']+'/DMS.example/3U5D.pdb')
    sasa_25S = format_data(pdb1, '1', Seq_25S.replace('\n',''))
    update_dms(sasa_25S, para['DataPath'], '25S.yeast')

if __name__ == '__main__':
    pymol.finish_launching()
    main_fun(main)
    pymol.cmd.quit()
