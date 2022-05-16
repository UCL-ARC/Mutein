"""
RSA 11.4.22

Script to take a file in the format given to me by Michael Hall and then produce a list of proteins

We use the bioservices python library to access databases
https://pypi.org/project/bioservices/

"""
import os
import sys
import yaml
import pandas as pd

dirs = os.path.dirname(os.path.realpath(__file__)).split("/")[:-2]
retpath = "/".join(dirs) + "/libs"
sys.path.append(retpath)
import Paths
import Arguments
import BatchMaker
import Gene
import Variant
import genetoprotein
import genestovariants
import SwissModel
import UniProt
import PdbRunner


def run_pipeline(args):
    argus = Arguments.Arguments(args)
    install_dir = argus.arg("install_dir")
    sys.path.append(install_dir)
    sys.path.append(install_dir + "/Pipelines")
    sys.path.append(install_dir + "/Pipelines/libs")
    data_dir = argus.arg("data_dir")
    dataset = argus.arg("dataset")
    dataset_path = Paths.Paths(data_dir, install_dir + "Pipelines/geneanalysis", dataset=dataset)
    genes = []
    """
    "gene_name"	"n_syn"	"n_mis"	"n_non"	"n_spl"	"n_ind"	"wmis_cv"	"wnon_cv"	"wspl_cv"	"wind_cv"	"pmis_cv"	"ptrunc_cv"	"pallsubs_cv"	"pind_cv"	"qmis_cv"	"qtrunc_cv"	"qallsubs_cv"	"pglobal_cv"	"qglobal_cv"
    "NOTCH1"	404	2700	772	360	1208	3.98846230912502	21.7266421919647	21.7266421919647	17.6547328391658	0	0	0	8.07429581976078e-68	0	0	0	0	0
    """
    # 1.) First prepare the variants file
    genes_variant_file = dataset_path.dataset_outputs + "genes_variants.txt"
    gene_variant_dic = genestovariants.extractVariantsFromFile(genes_variant_file)
    # 2.) Now go and find all the pdbs, also look at the variants per gene
    genes_file = dataset_path.dataset_outputs + "/genes.txt"
    batches = [] #all the gene batches for the entire dataset
    batches_stitch= [] #all the gene stitches for the entire dataset
    with open(genes_file, "r") as fr:
        lines = fr.readlines()
        #lines = ["","FAT1"]
        for l in range(1, len(lines)):
            # for l in range(1,2):
            line = lines[l]
            cols = line.split("\t")
            gene = cols[0].upper()            
            if gene[0] == '"':
                gene = gene[1:]
            if gene[-1] == '"':
                gene = gene[:-1]
            gene_path = Paths.Paths(data_dir,install_dir + "Pipelines/geneanalysis",dataset=dataset,gene=gene)
            accession = genetoprotein.accession_from_bioservices(gene)
            if len(accession) > 1:
                seq = genetoprotein.sequence_from_bioservices(accession)
                seq_lines = seq.split("\n")
                wholeseq = ""
                for s in range(1, len(seq_lines)):
                    sl = str(seq_lines[s].strip())
                    wholeseq += sl
                gn = Gene.Gene(gene, accession, wholeseq)
                genes.append(gn)  # main repository for data we are creating in function
                script_file = "libs/pipeline_qsubber.py"
                yaml_file = "geneanalysis/config/batch_pdb.yml"
                bm = BatchMaker.BatchMaker(script_file, yaml_file)
                # CREATE the variants for the gene
                vrs = gene_variant_dic[gene.upper()]
                for i in range(len(vrs["bases"])):
                    bs = vrs["bases"][i]
                    p1 = vrs["prev_aa"][i]
                    rs = vrs["residue"][i]
                    na = vrs["new_aa"][i]
                    vr = vrs["variant"][i]
                    # print(bs,p1,rs,na,vr)
                    vrnt = Variant.Variant(gene, vr, bs)
                    gn.addVariant(vrnt)
                # CREATE the pdbs for the gene                
                # both these searches retun a tuple list of the pdb code and the thruput gene file path, ready for pdb inputs
                up = UniProt.UniProt(gene)
                df, pdb_list_up = up.searchForStructures(gene_path.gene_outputs,gene_path.gene_outpdbs,fragment=50)
                accession = up.accession
                sm = SwissModel.SwissModel(gene,accession)
                dfs, pdb_list_sm = sm.searchForStructures(gene_path.gene_outputs,gene_path.gene_outpdbs,fragment=50)                
                dfs.append(df)
                pdb_list = pdb_list_up + pdb_list_sm
                vc = pd.concat(dfs, axis=0)
                vc.to_csv(gene_path.gene_outputs + "Coverage_all.csv",index=False)
                print(vc)                
                for pdb in pdb_list:                    
                    pdb_path = Paths.Paths(data_dir,install_dir + "Pipelines/geneanalysis",dataset=dataset,gene=gn.gene,pdb=pdb.pdbcode)
                    prun = PdbRunner.PdbRunner(pdb.pdbcode)
                    prun.copyToInput(gene_path.gene_outpdbs, pdb_path.pdb_inputs,vc)
                    dfp = gn.getPdbVariantCoverageDataFrame(pdb)
                    dfp.to_csv(pdb_path.pdb_inputs + "/variants.csv", index=False)
                    bm.addBatch(dataset, gn.gene, pdb.pdbcode)
                    gn.addPdb(pdb)
                                    
                # and we want only 1 batch for the stitching
                script_file = "libs/pipeline_qsubber.py"
                yaml_file = "geneanalysis/config/batch_genestitch.yml"
                bm2 = BatchMaker.BatchMaker(script_file, yaml_file)
                bm2.addBatch(dataset, gene, "x")
                bm2.printBatchScript(
                    gene_path.gene_outputs
                    + "/ppl_"
                    + dataset
                    + "_"
                    + gene
                    + "_stitch.sh",
                    "",
                )
                batches_stitch.append(bm2.printBatchScript(
                    gene_path.pipeline_path
                    + "/ppl_"
                    + dataset
                    + "_"
                    + gene
                    + "_stitch.sh",
                    gene_path.pipeline_path
                    + "/ppl_"
                    + dataset
                    + "_"
                    + gene
                    + "_STITCH_sym.sh",
                ))

            # having found our collection of genes with assopciated pdbs and variants we can now create the pdb datasets
            # for gn in genes:
            # gene_path = Paths.Paths("geneprot",dataset=dataset,gene=gn.gene)
            dfv = gn.getVariantCandidatesDataFrame()
            dfv.to_csv(gene_path.gene_outputs + "/pdb_candidates.csv", index=False)
            
            # And save it all in the dataset output
            # we do this each time so that if it gets abandoned we knw where we were
            dfp = genes[0].getDatasetGenesPdbsDataFrame(dataset, genes)
            dfp.to_csv(dataset_path.dataset_outputs + "/pdb_coverage.csv", index=False)

            bm.printBatchScript(
                gene_path.gene_outputs + "/ppl_" + dataset + "_" + gene + ".sh", ""
            )
            batch_to_add = bm.printBatchScript(
                gene_path.pipeline_path + "/ppl_" + dataset + "_" + gene + ".sh",
                gene_path.pipeline_path + "/ppl_" + dataset + "_" + gene + "sym.sh",
            )

            batches.append(batch_to_add)
        

    bmAll = BatchMaker.BatchMaker("", "")        
    bmAll.printBatches(
        gene_path.pipeline_path
        + "/ppl_"
        + dataset
        + "_all.sh",            
        batches)

    bmAll_stitch = BatchMaker.BatchMaker("", "")        
    bmAll_stitch.printBatches(
        gene_path.pipeline_path
        + "/ppl_"
        + dataset
        + "_all_stitch.sh",            
        batches_stitch)




    print("### COMPLETED genes to proteins pipeline ###")
    print("MUTEIN SCRIPT ENDED")


##########################################################################################
if __name__ == "__main__":
    import sys

    globals()["run_pipeline"](sys.argv)
