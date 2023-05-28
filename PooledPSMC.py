#!/usr/env/bin python
#####
import argparse

######
parser = argparse.ArgumentParser(description='PooledPSMC -- Auto-PSMC with bamlist as input and psmc results as output. This method pools all individuals and consider them as a same individual.')
parser.add_argument('--work_dir',help='working directory')
parser.add_argument('--name_list',help='a file that have all the project names (you can use individual names as project names). one each line')
parser.add_argument('--bam_list_file',help='a file that have all the bam path to individuals. one each line')
parser.add_argument('--ref_gen',help='reference genome path')
parser.add_argument('--d',default=6, help='only sites with depth over this are used')
parser.add_argument('--q',default=20, help='sequence quality; only sequence with quality over this is used')
parser.add_argument('--N',default=25, help='maximum number of iterations [25]')
parser.add_argument('--t',default=15, help='maximum 2N0 coalescent time [15]')
parser.add_argument('--r',default=5, help='initial theta/rho ratio [5]')
parser.add_argument('--p',default='64*1', help='pattern of parameters [64*1]')
parser.add_argument('--boostrap',default=100, help='default 100')
parser.add_argument('--T',default=1, help='thread; default 1')
parser.add_argument('--mu', help='mutation rate; for example: 1e-08')
parser.add_argument('--g', help='generation time (year); for example: 10')



args = parser.parse_args()
#######
import os
import subprocess
import time
from multiprocessing import Process, current_process, Lock, Pool
import sys
import pandas as pd
import matplotlib.pyplot as plt


########################################################
#### change dir
os.chdir(args.work_dir)

with open(f'../{args.bam_list_file}','r') as f:
    bam_list = [i.strip() for i in f.readlines() if not i.strip()=='']

with open(f'../{args.name_list}','r') as f:
    project_name_list = [i.strip() for i in f.readlines() if not i.strip()=='']

    
output_path=args.work_dir
reference_genome_path=args.ref_gen
filtering_depth = int(args.d)
sequence_quality = int(args.q)
N=int(args.N)   #maximum number of iterations [30]
t=int(args.t)    #maximum 2N0 coalescent time [15]
r=int(args.r)     #initial theta/rho ratio [4]
p=str(args.p)    #pattern of parameters [4+5*3+4]
boostrap=int(args.boostrap) #
thread=int(args.T)
mu=str(args.mu)
g=str(args.g)


########################################################
### define linux command wrapper
def process_command(command):
    p = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    while True:
        buff = p.stdout.readline()
        err = p.stderr.readline()
        buff = bytes.decode(buff)
        err = bytes.decode(err)
        print(buff,err)
        if buff == '' and p.poll() != None:
            break

    p.wait()
    
def path_check(path_list):
    for prereq in path_list:
        if not os.path.exists(prereq):
            print(f'{prereq} not exists')
            raise
            
            
def run_one_sample_until_psmcfa(param_dict, lock):
    
    #### read params
    project_name = param_dict['project_name']
    bam_path = param_dict['bam_path']
    output_path=param_dict['output_path']
    reference_genome_path = param_dict['reference_genome_path']
    filtering_depth = param_dict['filtering_depth']
    sequence_quality = param_dict['sequence_quality']
    N=param_dict['N']
    t=param_dict['t']
    r=param_dict['r']
    p=param_dict['p']
    boostrap=param_dict['boostrap']


    ### 1. combine sorted bam file
    command = f'samtools sort {bam_path}|bcftools mpileup -Ou -I -f {reference_genome_path} - |bcftools call -c -Ov -o {output_path}/{project_name}.vcf'
    with lock:
        print('Project name: ',project_name, '\nppid: ',os.getppid(), '\npid: ',os.getpid(),'\nCommand: ',command,'\n')
    process_command(command)
    with lock:
        print('done','\n')

    ### 2. make fastaq
    command = f'vcfutils.pl vcf2fq -d {filtering_depth} {output_path}/{project_name}.vcf| gzip > {output_path}/{project_name}.vcf.diploid.d{filtering_depth}.fq.gz'
    with lock:
        print('Project name: ',project_name, '\nppid',os.getppid(), '\npid',os.getpid(),'\nCommand: ',command,'\n')
    process_command(command)
    with lock:
        print('done','\n')

    ### 3. fq2psmcfa
    command = f'fq2psmcfa -q {sequence_quality} {output_path}/{project_name}.vcf.diploid.d{filtering_depth}.fq.gz > {output_path}/{project_name}.vcf.diploid.d{filtering_depth}.psmcfa'
    with lock:
        print('Project name: ',project_name, '\nppid',os.getppid(), '\npid',os.getpid(),'\nCommand: ',command,'\n')
    process_command(command)
    with lock:
        print('done','\n')



################## process each file
lock = Lock()
process_list=[]
for bam_path,project_name in zip(bam_list,project_name_list):
    params_dict={
        'project_name':project_name,
        'bam_path':bam_path,
        'output_path':output_path,
        'reference_genome_path':reference_genome_path,
        'filtering_depth':filtering_depth,
        'sequence_quality':sequence_quality,
        'N':N,   #maximum number of iterations [30]
        't':t,   #maximum 2N0 coalescent time [15]
        'r':r,    #initial theta/rho ratio [4]
        'p':p,     #pattern of parameters [4+5*3+4]
        'boostrap':boostrap #
    }
    
    ### append process
    process_list.append(Process(target=run_one_sample_until_psmcfa, args=(params_dict,lock)))
    
#### start processing
for process in process_list:
    process.start()
    
#### finish processing
for process in process_list:
    process.join()


### cat psmcfa (different chromosomes)
command = f'rm {output_path}/combined.psmcfa'
process_command(command)
for file in [i for i in os.listdir(output_path) if i.endswith('.psmcfa') and not i.startswith('combined')]:
    print(file)
    command = f'cat {output_path}/{file} >> {output_path}/combined.psmcfa'
    process_command(command)
    print('psmcfa concatenated')


### 4. run psmc
lock = Lock()
command = f'psmc -N{N} -t{t} -r{r} -p "{p}" -o {output_path}/combined.psmc {output_path}/combined.psmcfa'
with lock:
    print('\nppid',os.getppid(), '\npid',os.getpid(),'\nCommand: ',command,'\n')
process_command(command)
with lock:
    print('done','\n')
    
    
#### splitpsmcfa
lock = Lock()
command = f'splitfa {output_path}/combined.psmcfa > {output_path}/split.combined.psmcfa'

with lock:
    print('\nppid',os.getppid(), '\npid',os.getpid(),'\nCommand: ',command,'\n')
process_command(command)
with lock:
    print('done','\n')
    
            
#### boostraps
lock=Lock()

def run_a_single_boostrap(round_):
    command = f'psmc -N{N} -t{t} -r{r} -b -p "{p}" -o {output_path}/round{round_}split.combined.psmc {output_path}/split.combined.psmcfa'
    process_command(command)
    print(f'round {round_}: done')
    
with Pool(thread) as pool:
    # prepare arguments
    items = [(i,) for i in range(boostrap)]
    # issue tasks to the process pool and wait for tasks to complete
    pool.starmap(run_a_single_boostrap, items)
      
        
### cat and plot
command = f'psmc_plot.pl -R -u {mu} -g {g} -M "combined, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100" -p combined.boosts.mu{mu} {output_path}/combined.psmc {output_path}/round*split.combined.psmc'
# cat {output_path}/combined.psmc {output_path}/round*split.combined.psmc > combined.boosts.psmc 
lock=Lock()
with lock:
    print('\nppid: ',os.getppid(), '\npid: ',os.getpid(),'\nCommand: ',command,'\n')
process_command(command)
with lock:
    print('done','\n')
    


######### final plot
for file in [i for i in os.listdir('.') if i.endswith('.txt') and i.startswith('combined.boosts')]:
    data = pd.read_csv(file, sep='\t',header=None)
    data = data[data.iloc[:,0]>1e4]
    year_ = data[0].values
    ne_ = data[1].values
    plt.step(year_, ne_, linewidth=1, alpha=0.1, c='green')
    
    if file.split('.')[-2]=='0':
        plt.step(year_, ne_, linewidth=2, alpha=1, c='green')

        
plt.xscale("log")
plt.ylabel('Effective Population Size($10^{4}$)')
plt.xlabel('Time (years ago) g=10, mu=$10^{-8}$')
plt.tight_layout()
plt.savefig('Pooled_PSMC_results.pdf')
plt.show()


            
            
            
            
    

