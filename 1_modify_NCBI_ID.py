#### 1. 转换NCBI Database 中基因组数据的蛋白序列的 ID 以用于进化分析
### >lcl|CM022752.1_cds_KAF4093247.1_1 [protein=hypothetical protein] [protein_id=KAF4093247.1] [location=complement(join(8916..9148,9253..9471,9558..9902,10016..10357,13313..13582,13687..13905,13991..14335,14523..14577))] [gbkey=CDS]
###           to
### >KAF4093247.1_1
### 
### 2. 去除序列中多余的换行符
###  Usage: python3 1_modify_NCBI_ID.py FastaFile.fa


from os import write
import sys
import os


file = sys.argv[1]
path = os.path.abspath(file)
name =  path.replace(path.split("/")[-1], path.split("/")[-1].split('.')[0] + '.idfixed.fa')

print(f"Start converting ID format:\n  {name} \n")

with open(path) as f:
    for line in f:
        if line.startswith('>'):
            id = line.split(' ')
            id = id[0].split('_cds_')
            id = '>'+id[1]


            #id = id[-2].replace('[', '>') + '_'+ id[-1].replace('\n', '').replace(']', '|') + id[0].replace('>', '')
            #print(id + '\n')
            with open(name, 'a') as f:
                f.write('\n' + id + '\n')
        else:
            line = line.replace('\n', '')
            with open(name, 'a') as f:
                f.write(line)


print(f"Format conversion completed {name} \n")
                
