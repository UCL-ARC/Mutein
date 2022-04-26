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
import Pdb
import Variant

# import from the shared library in Mutein/Pipelines/shared/lib
dirs = os.path.dirname(os.path.realpath(__file__)).split("/")[:-2]
retpath = "/".join(dirs) + "/shared/libs"
sys.path.append(retpath)
import Paths
import Arguments



def run_pipeline(args):
    argus = Arguments.Arguments(args)
    dataset = argus.arg("dataset")
    dataset_path = Paths.Paths("vcf", dataset=dataset)
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
    with open(genes_file, "r") as fr:
        lines = fr.readlines()
        # lines = ["","FAT1"]
        for l in range(1, len(lines)):
            # for l in range(1,2):
            line = lines[l]
            cols = line.split("\t")
            gene = cols[0].upper()
            if gene[0] == '"':
                gene = gene[1:]
            if gene[-1] == '"':
                gene = gene[:-1]
            gene_path = Paths.Paths("geneprot", dataset=dataset, gene=gene)
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
                pdbs_paths = genetoprotein.pdbs_from_accession_bioservices(accession)
                print("gene=", gene, "accession=", accession)  # ,"pdbs=",pdbs_paths)
                for pdb_path in pdbs_paths:
                    pdb = pdb_path["pdb"]
                    url = pdb_path["path"]
                    biopdb = genetoprotein.retrievePdbStructure(
                        url, pdb, gene_path.gene_outpdbs + pdb + ".pdb"
                    )
                    if biopdb != None:
                        method, res = (
                            biopdb.header["structure_method"],
                            biopdb.header["resolution"],
                        )
                        import Bio.PDB as bio

                        ppb = (
                            bio.CaPPBuilder()
                        )  # PPBuilder is C-N and CAPPBuilder is CA-CA
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
                                    print("!!!", resis, gene, pdb)
                                # print(start,seq_one)
                                has_match = True
                                pb = Pdb.Pdb(
                                    gene,
                                    pdb,
                                    chain,
                                    start + 1,
                                    start + len(seq_one),
                                    method,
                                    res,
                                )
                                gn.addPdb(pb)
                        if not has_match:
                            genetoprotein.removePdbStructure(
                                url, pdb, gene_path.gene_outpdbs + pdb + ".pdb"
                            )

            # having found our collection of genes with assopciated pdbs and variants we can now create the pdb datasets
            # for gn in genes:
            # gene_path = Paths.Paths("geneprot",dataset=dataset,gene=gn.gene)
            dfv = gn.getVariantCandidatesDataFrame()
            dfv.to_csv(gene_path.gene_outputs + "/pdb_candidates.csv", index=False)
            dfp = gn.getPdbCoverageDataFrame()
            dfp.to_csv(gene_path.gene_outputs + "/pdb_coverage.csv", index=False)

            for pdbcod, pdb in gn.pdbs.items():
                # print(pdb.pdbcode)
                # only use x-ray and alphafold:
                if (
                    True
                ):  # pdb.getMethod() == "x-ray" or pdb.getMethod() == "alphafold" or pdb.getMethod() == "e-m":
                    pdb_path = Paths.Paths(
                        "pdb", dataset=dataset, gene=gn.gene, pdb=pdb.pdbcode
                    )
                    dfp = gn.getPdbVariantCoverageDataFrame(pdb)
                    dfp.to_csv(pdb_path.pdb_inputs + "/variants.csv", index=False)
                    pdb.downloadPdb(pdb_path.pdb_inputs)
                    # I don't think I need the config yaml anymore with the new batch system
                    #gn.createPdbConfigYaml(pdb, pdb_path.pdb_inputs + "/config.yml")

            # And save it all in the dataset output
            # we do this each time so that if it gets abandoned we knw where we were
            dfp = genes[0].getDatasetGenesPdbsDataFrame(dataset, genes)
            dfp.to_csv(dataset_path.dataset_outputs + "/pdb_coverage.csv", index=False)

    ##### INPUTS #############################################
    # The inputs to this function are the pdbfile and the chain id (might optionally consider the positionscan mutation type)
    # print("### Pipeline: genes to proteins ###")
    # path = Paths.Paths()
    # argus = Arguments.Arguments(args)
    # repair_path = argus.arg("interim_path") + "repair" + argus.arg("repairs") + "/"
    # argus.params["repair_path"] = repair_path
    # hlp.goto_job_dir(argus.arg("repair_path"), args, argus.params, "_inputs01")
    ############################################

    print("### COMPLETED genes to proteins pipeline ###")


##########################################################################################
if __name__ == "__main__":
    import sys

    globals()["run_pipeline"](sys.argv)
