"""
RSA 11.4.22

Script to take a file in the format given to me by Michael Hall and then produce a list of proteins

We use the bioservices python library to access databases
https://pypi.org/project/bioservices/

"""
import os
import sys
import genetoprotein
import genestovariants
import yaml
import pandas as pd
import Gene

#import from the shared library in Mutein/Pipelines/shared/lib
dirs = os.path.dirname(os.path.realpath(__file__)).split("/")[:-2]
retpath = "/".join(dirs) + '/shared/lib'
sys.path.append(retpath)
import Paths


def run_pipeline(args):
    dataset = 'shearwater'
    dataset_path = Paths.Paths('vcf',dataset=dataset)
    genes = []
    """
    "gene_name"	"n_syn"	"n_mis"	"n_non"	"n_spl"	"n_ind"	"wmis_cv"	"wnon_cv"	"wspl_cv"	"wind_cv"	"pmis_cv"	"ptrunc_cv"	"pallsubs_cv"	"pind_cv"	"qmis_cv"	"qtrunc_cv"	"qallsubs_cv"	"pglobal_cv"	"qglobal_cv"
    "NOTCH1"	404	2700	772	360	1208	3.98846230912502	21.7266421919647	21.7266421919647	17.6547328391658	0	0	0	8.07429581976078e-68	0	0	0	0	0
    """
    genes_file = dataset_path.dataset_outputs + '/genes.txt'
    genes_variant_file = dataset_path.dataset_outputs + '/genes_variants.txt'
    gene_variant_dic = genestovariants.extractVariantsFromFile(genes_variant_file)
    genes_paths_yaml = dataset_path.dataset_outputs + '/genes_pdbs.yml' # We want a yaml file for the whole dataset
    genes_variants_csv = dataset_path.dataset_outputs + '/genes_variants.csv' # 
    dic_genes_paths = {}    
    segments_list = []
    print(genes_file)
    with open(genes_file, 'r') as fr:                
        lines = fr.readlines()
        lines = ["notch1"]
        for line in lines:
            cols = line.split('\t')
            gene=cols[0]
            gene_path = Paths.Paths('geneprot',dataset=dataset,gene=gene)
            accession = genetoprotein.accession_from_bioservices(gene)            
            if len(accession)>1:
                seq = genetoprotein.sequence_from_bioservices(accession)
                seq_lines=seq.split('\n')
                wholeseq = ''
                for s in range(1,len(seq_lines)):
                    sl =  str(seq_lines[s].strip())
                    wholeseq += sl                    
                print(len(wholeseq),wholeseq)
                pdbs_paths = genetoprotein.pdbs_from_accession_bioservices(accession)
                dic_genes_paths[gene] = pdbs_paths
                pdb_paths_csv = gene_path.gene_outputs + '/pdbs_segments.csv' #we want an input yaml per gene
                print("gene=",gene, "accession=",accession)#,"pdbs=",pdbs_paths)                                                    
                pdbs_good = {}                              
                pdbs_good['gene'] = []
                pdbs_good['chain'] = []
                pdbs_good['segment'] = []                
                pdbs_good['pdb'] = []  
                pdbs_good['method'] = []  
                pdbs_good['resolution'] = []  
                #pdbs_good['url'] = []

                for pdb_path in pdbs_paths:
                    pdb = pdb_path['pdb']
                    url = pdb_path['path']
                    print(pdb,url)
                    biopdb = genetoprotein.retrievePdbStructure(url,pdb,gene_path.gene_outpdbs + "/" + pdb + ".pdb")
                    method,res = biopdb.header["structure_method"],biopdb.header["resolution"]
                    import Bio.PDB as bio
                    ppb = bio.CaPPBuilder() #PPBuilder is C-N and CAPPBuilder is CA-CA
                    has_match = False                    
                    for pp in ppb.build_peptides(biopdb):
                        seq_one = str(pp.get_sequence())                        
                        start = wholeseq.find(seq_one)
                        chain = ""
                        if start > -1:
                            try:
                                resis = pp.get_ca_list()[0]
                                chain = resis.parent.get_parent().id
                            except:                                
                                print("!!!",resis,gene,pdb)
                            print(start,seq_one)
                            has_match = True                            
                            segments_list.append([start,start+len(seq_one),pdb])
                            pdbs_good['gene'].append(gene)                            
                            pdbs_good['pdb'].append(pdb)                            
                            pdbs_good['segment'].append(str(start+1) + ":" + str(start+len(seq_one)))
                            pdbs_good['chain'].append(chain)
                            pdbs_good['method'].append(method)
                            pdbs_good['resolution'].append(res)
                            #pdbs_good['url'].append(url)
                    if not has_match:                                            
                        genetoprotein.removePdbStructure(url,pdb,gene_path.gene_outpdbs + "/" + pdb + ".pdb")
                
                pdbs_df = pd.DataFrame.from_dict(pdbs_good)
                pdbs_df.to_csv(pdb_paths_csv,index=False)                                
            else:
                print("gene=",gene, "accession=",accession)
    with open(genes_paths_yaml, 'w') as wf:
        outputs = yaml.dump(dic_genes_paths, wf)
    # now we know what the coverage is, we could go through the list of variants and tag which pdb files could be appropriate
    gene_variant_dic['candidates'] = []
    for rid in gene_variant_dic['residue']:
        candidate = []
        for start,end,pdb in segments_list:
            if int(start) <= int(rid) and int(end) >= int(rid):
                if pdb not in candidate:
                    candidate.append(pdb)
        gene_variant_dic['candidates'].append(" ".join(candidate))
    
    gene_variant_df = pd.DataFrame.from_dict(gene_variant_dic)
    gene_variant_df.to_csv(genes_variants_csv,index=False)
    
    #also save this in each gene folder
    genes = gene_variant_df["gene"].unique()
    for gene in genes:
        gene_one_df = gene_variant_df[gene_variant_df["gene"]==gene]
        gene_one_df = gene_one_df.sort_values(by='residue', ascending=True)
        gene_path = Paths.Paths('geneprot',dataset=dataset,gene=gene)        
        gene_one_df.to_csv(gene_path.gene_outputs + '/gene_variants.csv',index=False)

    # And now make a new protein folder ready for HPC and analysis
                # ....for the chosen pdb structures
                # ... ... which is currently only alphafold for the first pass
                for i in range(len(pdbs_df.index)):
                    pdb = pdbs_df['pdb'][i]
                    segment = pdbs_df['segment'][i]
                    start = int(segment.split(":")[0])
                    end = int(segment.split(":")[1])
                    chain = pdbs_df['chain'][i]

    








    
    
            

    
    
    
    
      

    ##### INPUTS #############################################
    # The inputs to this function are the pdbfile and the chain id (might optionally consider the positionscan mutation type)
    #print("### Pipeline: genes to proteins ###")
    #path = Paths.Paths()
    #argus = Arguments.Arguments(args)
    #repair_path = argus.arg("interim_path") + "repair" + argus.arg("repairs") + "/"
    #argus.params["repair_path"] = repair_path
    #hlp.goto_job_dir(argus.arg("repair_path"), args, argus.params, "_inputs01")
    ############################################

    print("### COMPLETED genes to proteins pipeline ###")


##########################################################################################
if __name__ == "__main__":
    import sys
    globals()["run_pipeline"](sys.argv)
