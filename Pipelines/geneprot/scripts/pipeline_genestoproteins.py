def run_pipeline(args):

    ##### INPUTS #############################################
    # The inputs to this function are the pdbfile and the chain id (might optionally consider the positionscan mutation type)
    print("### Pipeline: genes to proteins ###")
    argus = Arguments.Arguments(args)
    repair_path = argus.arg("interim_path") + "repair" + argus.arg("repairs") + "/"
    argus.params["repair_path"] = repair_path
    hlp.goto_job_dir(argus.arg("repair_path"), args, argus.params, "_inputs01")
    ############################################

    print("### COMPLETED genes to proteins pipeline ###")


##########################################################################################
if __name__ == "__main__":
    import sys

    globals()["run_pipeline"](sys.argv)
