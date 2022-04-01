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
import pandas as pd
from shutil import copyfile
import helper as hlp

def run_pipeline04(args):
    print('### Foldx aggregate pos scan ###')
    ##############################################
    iparams = hlp.inputparams(args)    
    pdb = ''
    if 'pdb' in iparams:
        pdb = iparams['pdb']
    cparams = hlp.configparams(pdb)    
    params = hlp.mergeparams(cparams,iparams)
    print(params)
    user = params['user']
    user, (foldxe, pythonexe, environment) = hlp.getenvironment(user)
    print(user, foldxe, pythonexe,environment)
    pdb = params['pdb']
    jobname = params['name']        
    input_path, thruput_path, interim_path, output_path = hlp.get_make_paths(pdb,jobname)
    ############################################    
    os.chdir(interim_path)

    params_file = interim_path + 'params.txt'
    rownum = 1
    with open(params_file) as fr:
        paramscontent = fr.readlines()        
        rownum = len(paramscontent)

    ddg_file = 'PS_' + pdb + '_rep_scanning_output.txt'
    with open(ddg_file,'w') as fw:
        for r in range(rownum):
            jobresults_file = interim_path + '/row' + str(r+1) + '/' + ddg_file
            print('row path',jobresults_file)
            if os.path.exists(jobresults_file):
                with open(jobresults_file) as fr:
                    jobcontent = fr.readlines()        
                fw.writelines(jobcontent)

    #Make a dataframe
    aa_dict = {'ALA':'A','CYS':'C','ASP':'D','GLU':'E','PHE':'F','GLY':'G','HIS':'H','ILE':'I','LYS':'K','LEU':'L','MET':'M','ASN':'N','PRO':'P','GLN':'Q','ARG':'R','SER':'S','THR':'T','VAL':'V','TRP':'W','TYR':'Y','H1S':'o','H2S':'e'}                 
    a_dict = {'A':'ALA','C':'CYS','D':'ASP','E':'GLU','F':'PHE','G':'GLY','H':'HIS','I':'ILE','K':'LYS','L':'LEU','M':'MET','N':'ASN','P':'PRO','Q':'GLN','R':'ARG','S':'SER','T':'THR','V':'VAL','W':'TRP','Y':'TYR','o':'H1S','e':'H2S'}
    ddg_dic = {'aa':[],'chain':[],'rid':[],'mut':[],'mutid':[],'ddg':[]}
    with open(ddg_file,'r') as fr:
        lines = fr.readlines()        
    for line in lines:
        if len(line) > 7:        
            lns = line.strip().split('\t')    
            print(lns)
            res = lns[0]    
            aa = res[:3]
            chain = res[3:4]
            rid = res[4:-1]
            mut = res[-1:]
            ddg = lns[1]
            print(aa,chain,rid,mut,ddg)
            ddg_dic['aa'].append(aa)
            ddg_dic['chain'].append(chain)
            ddg_dic['rid'].append(int(rid))
            ddg_dic['mut'].append(mut)
            ddg_dic['mutid'].append(aa_dict[aa]+chain + rid + mut)
            ddg_dic['ddg'].append(float(ddg))
    import pandas as pd
    ddg_df = pd.DataFrame.from_dict(ddg_dic)
    df_file = output_path + pdb + '_ddg_dataframe.csv'
    ddg_df.to_csv(df_file,index=False)

    #And save something visual as a starting point for some analysis
    import matplotlib.pyplot as plt
    import seaborn as sns
    fig,(ax1,ax2) = plt.subplots(1,2)
    fig.suptitle(pdb + ' background mutations\nddg <1=stabilising >2.5=destabilising')


    xax = 'rid'
    yax = 'mut'
    hue = 'ddg'

    # first plt
    sns.histplot(data=ddg_df,x='ddg',palette="tab20",ax=ax1,bins=50)
    ax1.set_ylabel('')
    #ax1.set_title('ddg >2.5 destabilising')


    vmin = -2.5#ddg_df[hue].min() #they have defined >2.5 as destabilising
    vmax = -1*vmin
    ddg_df = ddg_df.sort_values(by=yax, ascending=False)
    g = ax2.scatter(ddg_df[xax], ddg_df[yax], c=ddg_df[hue], cmap='Spectral', edgecolor='silver', alpha=0.35, linewidth=0.1, s=20,vmin=vmin,vmax=vmax)
    cb = plt.colorbar(g,extend='both')
    cb.set_label(hue)
    ax2.set_xlabel('residue no')
    ax2.set_ylabel('')
    #plt.legend(title=hue,bbox_to_anchor=(1.01, 1), loc=2, borderaxespad=0.,shadow=False,fancybox=False)  # Put the legend out of the figure

    plot_file = output_path + pdb + '_background_plot.png'
    plt.savefig(plot_file)
    print('### COMPLETED FoldX aggregate job ###')