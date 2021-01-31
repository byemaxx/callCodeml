#!/usr/bin/env python
# -*- coding: utf-8 -*-
# 1.获得最长密码子序列
# 2.在ID前加入物种名（取自文件名称）
import sys
import os

file = sys.argv[1]
path = os.path.abspath(file)
#path = ("/Users/byeamx/catfishes/pep_catfishes_kaks/Tachysurus_fulvidraco.faa")
name =  path.replace(path.split("/")[-1], path.split("/")[-1].split('.')[0] + '.fna'
)
species = path.split(r'/')[-1].replace('.faa', '')

print(f"开始精简 {species} \n")
re = {}  
with open(path) as f:
    for line in f:
        seq = []
        if line.startswith('>'):
            id = line.split('_')[:-1]
            id = str.join('_', id).replace('>','')
            id =  '>' + species + '|'+ id  
        else:
            seq.append(line)

        if id not in re:
            re[id] = seq
        else:
            re[id] += seq

maxseq = {}
for k,v in re.items():
    seq = max(v, key=len)
    maxseq[k] = seq

with open(name,'w') as f:
    for k,v in maxseq.items():
        f.write( k +'\n')
        tem = v.__str__()
        f.write(tem)
print(f"精简完成 {species} \n")
