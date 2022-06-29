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

import _helper
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
import FileDf


def run_pipeline(args):
    argus = Arguments.Arguments(args)
    install_dir = argus.arg("install_dir")    
    data_dir = argus.arg("data_dir")
    dataset = argus.arg("dataset", "")
    gene = argus.arg("gene")
    pdbs = 0
    for gene in [gene]:
        gene_path = Paths.Paths(
            data_dir, install_dir, dataset=dataset, gene=gene,readonly=False
        )
        # this contains the organism id from uniport, eg human=9606 and mouse=10090
        accession_path = gene_path.gene_inputs + "accessions.csv"
        acc_df = FileDf.FileDf(accession_path).openDataFrame()        
        idx0 = acc_df.index[0]
        organism_id = str(acc_df["organism"][idx0])
        accession = str(acc_df["accession"][idx0])
        
        if len(accession) > 1:            
            seq = genetoprotein.sequence_from_uniprot(accession)            
            gn = Gene.Gene(gene, accession, seq)
            # Read the variants from intray
            dfv = pd.read_csv(gene_path.gene_inputs + "/variants.csv")
            for i in range(len(dfv.index)):
                bs = dfv["bases"][i]
                vr = dfv["variant"][i]
                vrnt = Variant.Variant(gene, vr, bs)
                gn.addVariant(vrnt)
            # CREATE the pdbs for the gene
            # both these searches retun a tuple list of the pdb code and the thruput gene file path, ready for pdb inputs
            frag = -1
            up = UniProt.UniProt(gene.upper(),accession)
            df, pdb_list_up = up.searchForStructures(
                gene_path.gene_outputs, gene_path.gene_outpdbs, fragment=frag
            )
            accession = up.accession
            sm = SwissModel.SwissModel(gene, accession)
            dfs, pdb_list_sm = sm.searchForStructures(
                gene_path.gene_outputs, gene_path.gene_outpdbs, fragment=frag
            )
            dfs.append(df)
            pdb_list = pdb_list_up + pdb_list_sm
            vc = pd.concat(dfs, axis=0)
            vc.to_csv(gene_path.gene_outputs + "Coverage_all.csv", index=False)
            print(vc)
            for pdb in pdb_list:
                pdb_path = Paths.Paths(
                    data_dir,
                    install_dir,
                    dataset=dataset,
                    gene=gn.gene,
                    pdb=pdb.pdbcode,
                    readonly=False
                )
                prun = PdbRunner.PdbRunner(pdb.pdbcode)
                prun.copyToInput(gene_path.gene_outpdbs, pdb_path.pdb_inputs, vc)
                dfp = gn.getPdbVariantCoverageDataFrame(pdb)
                dfp.to_csv(pdb_path.pdb_inputs + "/variants.csv", index=False)
                # bm.addBatch(dataset, gn.gene, pdb.pdbcode)
                gn.addPdb(pdb)
                pdbs += 1

            dfpdb = gn.getPdbTaskList()
            dfpdb.to_csv(gene_path.gene_outputs + "/pdb_tasklist.csv", index=False)

            # having found our collection of genes with assopciated pdbs and variants we can now create the pdb datasets
            # for gn in genes:
            # gene_path = Paths.Paths("geneprot",dataset=dataset,gene=gn.gene)
            dfv = gn.getVariantCandidatesDataFrame()
            dfv.to_csv(gene_path.gene_outputs + "/pdb_candidates.csv", index=False)

    print("### COMPLETED gene to proteins pipeline ###")
    #print("MUTEIN SCRIPT ENDED")
    return pdbs


##########################################################################################
if __name__ == "__main__":
    import sys

    globals()["run_pipeline"](sys.argv)
