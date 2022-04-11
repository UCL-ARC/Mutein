"""
RSA 11.4.22

Script to take a file in the format given to me by Michael Hall and then produce a list of proteins

We use the bioservices python library to access databases
https://pypi.org/project/bioservices/

"""
import os
import sys
import genetoprotein

#import from the shared library in Mutein/Pipelines/shared/lib
dirs = os.path.dirname(os.path.realpath(__file__)).split("/")[:-2]
retpath = "/".join(dirs) + '/shared/lib'
sys.path.append(retpath)
import Paths


def run_pipeline(args):
    path = Paths.Paths(['geneprot'],dataset='shearwater',gene='notch1')
    """
    "gene_name"	"n_syn"	"n_mis"	"n_non"	"n_spl"	"n_ind"	"wmis_cv"	"wnon_cv"	"wspl_cv"	"wind_cv"	"pmis_cv"	"ptrunc_cv"	"pallsubs_cv"	"pind_cv"	"qmis_cv"	"qtrunc_cv"	"qallsubs_cv"	"pglobal_cv"	"qglobal_cv"
    "NOTCH1"	404	2700	772	360	1208	3.98846230912502	21.7266421919647	21.7266421919647	17.6547328391658	0	0	0	8.07429581976078e-68	0	0	0	0	0
    """
    genes_file = path.dataset_outputs + '/genes.txt'
    genes_variant_file = path.dataset_outputs + '/genes_variants.txt'
    print(genes_file)
    with open(genes_file, 'r') as fr:                
        lines = fr.readlines()
        #lines = ["notch1"]
        for line in lines:
            cols = line.split('\t')
            gene=cols[0]            
            accession = genetoprotein.accession_from_bioservices(gene)
            if len(accession)>1:
                pdbs = genetoprotein.pdbs_from_accession_bioservices(accession)
                print("gene=",gene, "accession=",accession,"pdbs=",pdbs)
            else:
                print("gene=",gene, "accession=",accession)
            

    
    
    
    
      

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
