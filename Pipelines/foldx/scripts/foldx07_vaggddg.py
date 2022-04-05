'''
-----------------------------
RSA 17/03/22
-----------------------------

This aggregates the outputs from split positionscans into 1 file
-----------------------------
N.b this file may be run on the myriad clusters or on a local machine
-----------------------------
'''
import os
import statistics
import pandas as pd
from shutil import copyfile
import helper as hlp

def run_pipeline07(args):
    ##### INPUTS #############################################
    print('### Foldx variant aggregate ###')
    ##############################################
    iparams = hlp.inputparams(args)    
    pdb = ''
    if 'pdb' in iparams:
        pdb = iparams['pdb']
    cparams = hlp.configparams(pdb)    
    params = hlp.mergeparams(cparams,iparams)
    print('FINAL PARAMS',params)
    user = params['user']
    user, (foldxe, pythonexe, environment) = hlp.getenvironment(user)
    print(user, foldxe, pythonexe,environment)
    pdb = params['pdb']
    repairs = params['repairs']
    reppdb = params['pdb'] + '_rep' + repairs
    jobname = params['name']        
    input_path, thruput_path, interim_path, output_path = hlp.get_make_paths(pdb,jobname)
    vagg_path = interim_path + 'vagg/'
    hlp.goto_job_dir(vagg_path,args,params,'_inputs07') 

    params_file = thruput_path + 'variant_params.txt'

    variant_dirs = []
    with open(params_file) as fr:
        paramscontent = fr.readlines()
        for pc in paramscontent:        
            pcs = pc.strip().split(' ')  
            vnam = pcs[2].replace(',','_')
            vdir = pcs[3]
            variant_dirs.append([vdir,vnam])            

    ddg_dic = {'mutid':[],'ddg':[],'tag':[]}
    for vd,vn in variant_dirs:
        print(vd,vn)
        ddg_file = interim_path + vd + '/Dif_' + reppdb + '.fxout' #the pdb repaired files are always pdbcode_rep                
        print(ddg_file)
        if os.path.exists(ddg_file):
            with open(ddg_file) as fr:
                jobcontent = fr.readlines()
                energy_list = []            
                for linecontents in jobcontent:
                    line = linecontents.split('\t')
                    if line[0].endswith('.pdb'):
                        energy = line[1]
                        energy_list.append(float(energy))
                DDG = statistics.mean(energy_list)
                ddg_dic['mutid'].append(vn)
                ddg_dic['ddg'].append(DDG)
                ddg_dic['tag'].append(vd)
                
    #Make a dataframe
    import pandas as pd
    ddg_df = pd.DataFrame.from_dict(ddg_dic)
    df_file = reppdb + '_variants_ddg_dataframe.csv'
    ddg_df.to_csv(output_path+df_file,index=False)
    print('saved dataframe to',output_path + df_file)

    #And save something visual as a starting point for some analysis
    import matplotlib.pyplot as plt
    import seaborn as sns    
    fig,(ax1,ax2) = plt.subplots(1,2,figsize=(8,9))
    fig.suptitle(pdb + ' variant mutations\nddg <-1=stabilising >2.5=destabilising')    
    sns.set_color_codes("pastel")
    # first plt
    sns.barplot(x="ddg", y="tag", data=ddg_df,color="g",ax=ax1)
    ax1.set_ylabel('')
    # second plt
    sns.histplot(data=ddg_df,x='ddg',palette="tab20",ax=ax2,bins=20)
    ax2.set_ylabel('')
    sns.despine(left=True, bottom=True)
    plt.rcParams["axes.labelsize"] = 25
    plot_file = reppdb + '_variant_plot.png'
    plt.savefig(output_path+plot_file)
    print('saved plot to',output_path+plot_file)
    
    print('### COMPLETED FoldX aggregate job ###')
#####################################################################################################
if __name__ == '__main__':
    import sys
    globals()['run_pipeline07'](sys.argv)