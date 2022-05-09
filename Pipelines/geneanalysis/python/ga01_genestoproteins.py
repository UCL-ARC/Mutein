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

# import from the shared library in Mutein/Pipelines/shared/lib
dirs = os.path.dirname(os.path.realpath(__file__)).split("/")[:-2]
retpath = "/".join(dirs) + "/shared/libs"
sys.path.append(retpath)
import Paths
import Arguments
import BatchMaker
import Gene
import Pdb
import Variant
import genetoprotein
import genestovariants


def run_pipeline(args):
    argus = Arguments.Arguments(args)
    install_dir = argus.arg("install_dir")
    data_dir = argus.arg("data_dir")
    dataset = argus.arg("dataset")
    dataset_path = Paths.Paths("vcf",data_dir, install_dir+"Pipelines/geneanalysis", dataset=dataset)
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
            gene_path = Paths.Paths("geneprot", data_dir,install_dir+"Pipelines/geneanalysis",dataset=dataset, gene=gene)
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
                script_file = "lib/pipeline_qsubber.py"
                yaml_file = "geneanalysis/config/batch_pdb.yml"
                bm = BatchMaker.BatchMaker(script_file,yaml_file)
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
                        method, reso = (
                            biopdb.header["structure_method"],
                            biopdb.header["resolution"],
                        )
                        res = 0
                        try:
                            res = float(reso)
                        except:
                            pass

                        import Bio.PDB as bio

                        ppb = (
                            bio.CaPPBuilder()
                        )  # PPBuilder is C-N and CAPPBuilder is CA-CA
                        has_match = False
                        peptides = []
                        try:
                            peptides = ppb.build_peptides(biopdb)
                        except:
                            print("Build peptides error", pdb)

                        for pp in peptides:
                            seq_one = str(pp.get_sequence())
                            start = wholeseq.find(seq_one)
                            start += 1
                            chain = ""
                            if start > 0:
                                if True:
                                    # try:
                                    resis = pp.get_ca_list()[0]
                                    chain = resis.parent.get_parent().id
                                    # https://biopython.org/docs/1.75/api/Bio.PDB.Atom.html
                                    # get_full_id(self): Return the full id of the atom.
                                    # The full id of an atom is the tuple (structure id, model id, chain id, residue id, atom name, altloc).
                                    residue_num = resis.get_full_id()[3][1]
                                    offset = start - residue_num
                                    print(
                                        pdb,
                                        "gene=",
                                        start,
                                        "pdb=",
                                        residue_num,
                                        "offset=",
                                        offset,
                                        seq_one,
                                    )
                                    has_match = True
                                    pb = gn.getPdb(pdb, method, res)
                                    pb.addSegment(
                                        start, start - 1 + len(seq_one), offset, chain
                                    )
                                # except:
                                #    print("!!!", resis, gene, pdb)
                        if not has_match:
                            genetoprotein.removePdbStructure(
                                url, pdb, gene_path.gene_outpdbs + pdb + ".pdb"
                            )
                                     
                                            
                # and we want only 1 batch for the stitching
                script_file = "lib/pipeline_qsubber.py"
                yaml_file = "geneanalysis/config/batch_genestitch.yml"
                bm2 = BatchMaker.BatchMaker(script_file,yaml_file)                
                bm2.addBatch(dataset, gene, "x")
                bm2.printBatchScript(gene_path.gene_outputs+"/ppl_"+dataset+"_"+gene+"_stitch.sh","")
                bm2.printBatchScript(gene_path.pipeline_path+"/ppl_"+dataset+"_"+gene+"_stitch.sh",gene_path.pipeline_path+"/ppl_"+dataset+"_"+gene+"_STITCH")

            # having found our collection of genes with assopciated pdbs and variants we can now create the pdb datasets
            # for gn in genes:
            # gene_path = Paths.Paths("geneprot",dataset=dataset,gene=gn.gene)
            dfv = gn.getVariantCandidatesDataFrame()
            dfv.to_csv(gene_path.gene_outputs + "/pdb_candidates.csv", index=False)
            dfp = gn.getPdbCoverageDataFrame()
            dfp.to_csv(gene_path.gene_outputs + "/pdb_coverage.csv", index=False)

            for pdbcod, pdb in gn.pdbs.items():
                print("Gene contains:", pdb.pdbcode)
                # only use x-ray and alphafold:
                if pdb.getMethod() != "nmr":
                    pdb_path = Paths.Paths("pdb", data_dir,install_dir+"Pipelines/geneanalysis",dataset=dataset, gene=gn.gene, pdb=pdb.pdbcode)                    
                    dfp = gn.getPdbVariantCoverageDataFrame(pdb)
                    dfp.to_csv(pdb_path.pdb_inputs + "/variants.csv", index=False)
                    pdb.downloadPdb(pdb_path.pdb_inputs)
                    dfp = gn.getSinglePdbCoverageDataFrame(pdb)
                    dfp.to_csv(pdb_path.pdb_inputs + "/coverage.csv", index=False)

                    # so we want to use it in analysis
                    bm.addBatch(dataset, gn.gene, pdbcod)

            # And save it all in the dataset output
            # we do this each time so that if it gets abandoned we knw where we were
            dfp = genes[0].getDatasetGenesPdbsDataFrame(dataset, genes)
            dfp.to_csv(dataset_path.dataset_outputs + "/pdb_coverage.csv", index=False)

            bm.printBatchScript(gene_path.gene_outputs+"/ppl_"+dataset+"_"+gene+".sh","")
            bm.printBatchScript(gene_path.pipeline_path +"/ppl_"+dataset+"_"+gene+".sh",gene_path.pipeline_path +"/ppl_"+dataset+"_"+gene)

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
