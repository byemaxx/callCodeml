#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# callCodeML 
# Verson: 0.2
# Date: 2020.12.14
# Author:  Wu Qing
# This is used to call Codeml to calculate the positive selection of gene in Branch-Sites Model.(Require codon sequence)
# Usage: python3 callCodeml.py seqDir treeFile

import re
import sys
import os
import subprocess
from scipy.stats import chi2
import time
import shutil
from tqdm import tqdm


def help_info():
    print("This script is used to call codeML.\nUsage: python3 callCodeml.py Dir treeFile\n")

def getArgs():
    global dir_path
    global seqs
    global tree_file
    global start_time
    start_time = time.strftime('%Y%m%d-%H%M%S') # to keep the working dircitory name
    args = sys.argv[1:]
    if len(args) == 0:
        help_info()
    if os.path.isdir(args[0]) is True:
        dir_path = os.path.abspath(args[0])
        seqs = os.listdir(dir_path)
        seqs = [x for x in seqs if not str(x).startswith('.') ]
        tree_file = os.path.abspath(args[1])
        print(f"The drictory is: {dir_path}\nThe treeFile is: {tree_file}")
        #os.chdir(os.path.abspath(args[0]))
        return (dir_path, seqs, tree_file)
    else:
        help_info()




path0 = os.getcwd()



def creat_ctl(path_seq, name):

    seqfile = path_seq
    treefile = tree_file
       
    ctl_null = f'''
    seqfile = {seqfile}
    treefile = {treefile}
    outfile = {path0}/WorkingDrictory_{start_time}/{name}/null/{name}_null.res

    noisy = 9   
    verbose = 0   
    runmode = 0  

    seqtype = 1   
    CodonFreq = 2   
    clock = 0   
    aaDist = 0
    model = 2

    NSsites = 2   
    icode = 0   
    Mgene = 0

    fix_kappa = 0   
    kappa = 2   
    fix_omega = 1   
    omega = 1   

    fix_alpha = 1   
    alpha = .0  
    Malpha = 0   
    ncatG = 3   

    getSE = 0   
    RateAncestor = 0   

    fix_blength = 0  
    method = 0   
    Small_Diff = .45e-6
    cleandata = 1
    '''
    ctl_alte = ctl_null.replace(f'{path0}/WorkingDrictory_{start_time}/{name}/null/{name}_null.res', 
    f'{path0}/WorkingDrictory_{start_time}/{name}/alte/{name}_alte.res').replace('fix_omega = 1', 'fix_omega = 0').replace('omega = 1', 'omega = 2')

    with open(f'{path0}/WorkingDrictory_{start_time}/{name}/null/null.ctl', 'w') as f:
        f.write(ctl_null)
    
    with open(f'{path0}/WorkingDrictory_{start_time}/{name}/alte/alte.ctl', 'w') as f:
        f.write(ctl_alte)

def create_dir(name):
    try:
        os.mkdir(f'{path0}/WorkingDrictory_{start_time}/')
    except:pass
    try:
        os.mkdir(f'{path0}/WorkingDrictory_{start_time}/{name}')
    except:pass
    try:
        os.mkdir(f'{path0}/WorkingDrictory_{start_time}/{name}/null')
        os.mkdir(f'{path0}/WorkingDrictory_{start_time}/{name}/alte')
    except:pass

def call_codeml(type, path_seq, name):

    os.chdir(f'{path0}/WorkingDrictory_{start_time}/{name}/{type}')
    print(f"\n{time.strftime('%Y-%m-%d %H:%M:%S')} Staring  ModelA {type} of {name} ...")
    b = subprocess.Popen(['codeml', f'{path0}/WorkingDrictory_{start_time}/{name}/{type}/{type}.ctl'], stdout=subprocess.PIPE)
    b.wait()
    if b.returncode == 0:
        with open(f"{path0}/WorkingDrictory_{start_time}/{name}/{type}/{name}_{type}.res", "r") as f:
            t =  f.read()
            kappa = re.findall(r"kappa.+", t)[0].split(' ')[-1].replace('\n', '')
            ls = re.findall(r'lnL.+', t)[0].split(' ')
            ls=[x for x in ls if x!='']
            lnL = ls[4]
            np = ls[3].replace("):","")

            os.chdir(path0)
            return([lnL, np, kappa])            
    else:
        print(f"{time.strftime('%Y-%m-%d %H:%M:%S')} ERRO in run  ModelA {type} of {path_seq}\n")
    
    
    

def run(path_seq, name):
    create_dir(name)
    creat_ctl(path_seq, name)

    res_null = call_codeml('null', path_seq, name)
    res_alte = call_codeml('alte', path_seq, name)

    lnL0 = float(res_null[0])
    lnL1 = float(res_alte[0])
    np0 = int(res_null[1])
    np1 = int(res_alte[1])
    kappa = float(res_alte[2])

    lnl = abs(lnL0 - lnL1 )*2
    np = abs(np0 - np1)
    pvalue = 1 - chi2.cdf(lnl, np)

    print(f"\n{time.strftime('%Y-%m-%d %H:%M:%S')}\nName: {name} \nlnL0: {lnL0} lnL1: {lnL1} \nnp0: {np0} np1: {np1} \nKappa: {kappa} P-Value: {pvalue}\n")

    if pvalue < 0.05:
        try:
            os.mkdir(f'{path0}/result_{start_time}')
        except:
            pass
        #shutil.copy(f'WorkingDrictory_{start_time}/{name}/null/{name}_null.res', f'result_{start_time}/')
        shutil.copy(f'{path0}/WorkingDrictory_{start_time}/{name}/alte/{name}_alte.res', f'result_{start_time}/')

        with open(f'{path0}/result_{start_time}/result.tsv', 'a') as f:
            f.write(f"{name}\t{lnL0}\t{lnL1}\t{np0}\t{np1}\t{kappa}\t{pvalue}\n")
        


    


def main():
    t0 = time.time()
    try:
        getArgs()
        for i in tqdm(seqs):
            path_seq = dir_path + '/' + i
            name =  i.split('.')[0]
            run(path_seq, name)

            # t = threading.Thread(target=run, args={path_seq, name})
            # t.start()
            
            # p = multiprocessing.Process(target=run, args={path_seq, name})
            # p.start()
    
    finally:
        pass
    
    print(f'Time used: {time.time() - t0} S\n')

main()



