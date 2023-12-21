


codons_aa = {}

#TT
codons_aa["TTT"]="F"
codons_aa["TTC"]="F"
codons_aa["TTA"]="L"
codons_aa["TTG"]="L"
#CT
codons_aa["CTT"]="L"
codons_aa["CTC"]="L"
codons_aa["CTA"]="L"
codons_aa["CTG"]="L"
#AT
codons_aa["ATT"]="I"
codons_aa["ATC"]="I"
codons_aa["ATA"]="I"
codons_aa["ATG"]="M"
#GT
codons_aa["GTT"]="V"
codons_aa["GTC"]="V"
codons_aa["GTA"]="V"
codons_aa["GTG"]="V"
#TC
codons_aa["TCT"]="S"
codons_aa["TCC"]="S"
codons_aa["TCA"]="S"
codons_aa["TCG"]="S"
#CC
codons_aa["CCT"]="P"
codons_aa["CCC"]="P"
codons_aa["CCA"]="P"
codons_aa["CCG"]="P"
#AC
codons_aa["ACT"]="T"
codons_aa["ACC"]="T"
codons_aa["ACA"]="T"
codons_aa["ACG"]="T"
#GC
codons_aa["GCT"]="A"
codons_aa["GCC"]="A"
codons_aa["GCA"]="A"
codons_aa["GCG"]="A"
#TA
codons_aa["TAT"]="Y"
codons_aa["TAC"]="Y"
codons_aa["TAA"]="X"
codons_aa["TAG"]="X"
#CA
codons_aa["CAT"]="H"
codons_aa["CAC"]="H"
codons_aa["CAA"]="Q"
codons_aa["CAG"]="Q"
#AA
codons_aa["AAT"]="N"
codons_aa["AAC"]="N"
codons_aa["AAA"]="K"
codons_aa["AAG"]="K"
#GA
codons_aa["GAT"]="D"
codons_aa["GAC"]="D"
codons_aa["GAA"]="E"
codons_aa["GAG"]="E"
#TG
codons_aa["TGT"]="C"
codons_aa["TGC"]="C"
codons_aa["TGA"]="X"
codons_aa["TGG"]="W"
#CG
codons_aa["CGT"]="R"
codons_aa["CGC"]="R"
codons_aa["CGA"]="R"
codons_aa["CGG"]="R"
#AG
codons_aa["AGT"]="S"
codons_aa["AGC"]="S"
codons_aa["AGA"]="R"
codons_aa["AGG"]="R"
#GG
codons_aa["GGT"]="G"
codons_aa["GGC"]="G"
codons_aa["GGA"]="G"
codons_aa["GGG"]="G"

def getAA(triple_seq,append_first="", phase=0):
    #print(triple_seq,append_first,phase)    
    seq_out = ""
    end = ""        
    seq_use = append_first + triple_seq[int(phase):]
    for i in range(0,len(seq_use),3):
        nuc = seq_use[i:i+3].upper()
        if nuc in codons_aa:
            seq_out += codons_aa[nuc]
        else:
            #seq_out += "?" #we pass the end back, don;t need to add it on.
            end = nuc
            #raise Exception("Not a multiple of 3")
    return seq_use,seq_out,end


    