# snakemake 


# Field variables
configfile: "config/config.yaml"

# Local plan
# n: 1 2 4 8 16 32 45
# t: 1 2 4 8 16 32 64 128
# 3 replicates: "a", "b", "c"



"""

# 32
snakemake --profile profile/default --config n=2 t=32 s=methanobrevibacteraceae
snakemake --profile profile/default --config n=4 t=32 s=methanobrevibacteraceae
snakemake --profile profile/default --config n=8 t=32 s=methanobrevibacteraceae
snakemake --profile profile/default --config n=16 t=32 s=methanobrevibacteraceae
snakemake --profile profile/default --config n=32 t=32 s=methanobrevibacteraceae
snakemake --profile profile/default --config n=44 t=32 s=methanobrevibacteraceae

# 32
snakemake --profile profile/default --config n=2 t=32 s=prevotella
snakemake --profile profile/default --config n=4 t=32 s=prevotella
snakemake --profile profile/default --config n=8 t=32 s=prevotella
snakemake --profile profile/default --config n=16 t=32 s=prevotella
snakemake --profile profile/default --config n=32 t=32 s=prevotella
snakemake --profile profile/default --config n=64 t=32 s=prevotella
snakemake --profile profile/default --config n=124 t=32 s=prevotella


#  64 (Maybe not realistic to run with 64 threads as others are also using the machine.)


"""

# HPC plan
# I'm not sure it makes any sense to run this on a HPC as the other tools are already slow and running on a HPC would just make them look even worse. Maybe I can mention that in the paper.

# 0.00user        0.00system      0:00.10elapsed  0%CPU   (0text+0data    1952max)k       0inputs+0outputs        (0major+95minor)pagefaults      0swaps
# /usr/bin/time -f "%U\t%S\t%E\t%P\t%X\t%D\t%M\t%I\t%O\t%F\t%R\t%W" sleep 0.1
time_format = "-f \"%U\\t%S\\t%E\\t%P\\t%X\\t%D\\t%M\\t%I\\t%O\\t%F\\t%R\\t%W\""
columns = "user system elapsed cpu text_k data_k max_k inputs outputs major_pf minor_pf swaps"

# Defining temp as a function that does nothing, is a neat way of disabling deletion of output dirs.
#temp = lambda x: x

print(f"n{config['n']}_t{config['t']}_s{config['s']}")


# Sorry for mixing genus and family on this one.
source = config["s"] # prevotella | methanobrevibacteraceae



replicates = ["a", "b", "c"] # How many is enough, it is just to have some indication that the running time is stable.
#replicates = ["a"] # QUICK/DEBUG

rule all:
    input: 
        expand("results/time_tormes_n{n}_t{t}_r{r}_s{s}.txt", n = config['n'], t = config['t'], r = replicates, s = source),
        expand("results/time_asscom2_n{n}_t{t}_r{r}_s{s}.txt", n = config['n'], t = config['t'], r = replicates, s = source),
        expand("results/time_bactopia_n{n}_t{t}_r{r}_s{s}.txt", n = config['n'], t = config['t'], r = replicates, s = source),
        expand("results/time_asscom2bt_n{n}_t{t}_r{r}_s{s}.txt", n = config['n'], t = config['t'], r = replicates, s = source),


rule tormes:
    output:
        flag = touch("results/time_tormes_n{n}_t{t}_r{r}_s{s}.txt"),
        dir = directory("res_tormes_n{n}_t{t}_r{r}_s{s}"),
    benchmark: "benchmarks/tormes_n{n}_t{t}_r{r}_s{s}.tsv"
    conda: "envs/tormes.yaml"
    threads: 1 # This is just to force it to spawn only one process at a time
    shell: """
        
        # Only once run.
        # tormes-setup
        
        # Produce metadata with correct number of samples.
        # We can't write to the output dir yet wtf so we have to move that somewhere else wtf.
        # Plus one for the header.
        head -n $(( {wildcards.n}+1 )) assets/tormes_metadata_{wildcards.s}.tsv > .tmp_tormes_metadata_n{wildcards.n}.tsv
        
        # Tormes expects an empty directory, so we must first delete the one that snakemake makes ready.
        test -d {output.dir} && rm -r {output.dir}

        tormes \
            --metadata .tmp_tormes_metadata_n{wildcards.n}.tsv \
            --threads {wildcards.t} \
            --output {output.dir}
            
    """


# 1 Run asscom2 in "tormes" mode
rule asscom2:
    output:
        flag = touch("results/time_asscom2_n{n}_t{t}_r{r}_s{s}.txt"),
        dir = temp(directory("res_asscom2_n{n}_t{t}_r{r}_s{s}")),
    benchmark: "benchmarks/asscom2_n{n}_t{t}_r{r}_s{s}.tsv"
    conda: "envs/asscom2.yaml"
    threads: 1 # This is just to force it to spawn only one process at a time. Should be set as standard for all in the profiles instead.
    shell: """
       
        
        mkdir -p {output.dir}
        echo 1
        
        # Create fofn
        realpath {wildcards.s}_fna/*.fna > temp.txt
        head -n {wildcards.n} temp.txt > {output.dir}/fofn.txt
        
        echo 2

        head {output.dir}/fofn.txt

        
        export ASSCOM2_BASE="$(realpath ~/asscom2)"
        export ASSCOM2_PROFILE="${{ASSCOM2_BASE}}/profile/conda/default"
        
        
        ${{ASSCOM2_BASE}}/asscom2 \
            --cores {wildcards.t} \
            --config \
                fofn={output.dir}/fofn.txt \
                prokka_rfam=false \
                output_directory={output.dir} \
                prokka_kingdom=archaea \
            --until prokka abricate assembly_stats mlst panaroo gtdbtk
    

    """
    
    
    
rule bactopia:
    output:
        flag = touch("results/time_bactopia_n{n}_t{t}_r{r}_s{s}.txt"),
        dir = temp(directory("res_bactopia_n{n}_t{t}_r{r}_s{s}")),
    #conda: "envs/bactopia.yaml"
    benchmark: "benchmarks/bactopia_n{n}_t{t}_r{r}_s{s}.tsv"
    threads: 1 # This is just to force it to spawn only one process at a time
    shell: """
        
        source activate /glittertind/home/carl/miniforge3/envs/bactopia
        
        # Produce metadata table using built in "prepare" functionality
        # Plus one for the header.
        bactopia prepare --path {wildcards.s}_fna --assembly-ext .fna > .tmp_bactopia_fofn_nall.txt
        
        head -n $(( {wildcards.n}+1 )) .tmp_bactopia_fofn_nall.txt > .tmp_bactopia_fofn_n{wildcards.n}.txt

        bactopia  \
            --samples .tmp_bactopia_fofn_n{wildcards.n}.txt \
            --max_cpus 4 \
            -qs $(({wildcards.t} / 4))
            
        # Move results out of the way.
        mv bactopia {output.dir}
        mv work {output.dir}

    """


# 2 Run asscom2 in "bactopia" mode
rule asscom2bt:
    output:
        flag = touch("results/time_asscom2bt_n{n}_t{t}_r{r}_s{s}.txt"),
        dir = temp(directory("res_asscom2bt_n{n}_t{t}_r{r}_s{s}")),
    conda: "envs/asscom2.yaml"
    benchmark: "benchmarks/asscom2bt_n{n}_t{t}_r{r}_s{s}.tsv"
    threads: 1 # This is just to force it to spawn only one process at a time
    shell: """
        
        mkdir -p {output.dir}
        echo 1
        
        # Create fofn
        realpath {wildcards.s}_fna/*.fna > temp.txt
        head -n {wildcards.n} temp.txt > {output.dir}/fofn.txt
        
        echo 2

        head {output.dir}/fofn.txt

        
        export ASSCOM2_BASE="$(realpath ~/asscom2)"
        export ASSCOM2_PROFILE="${{ASSCOM2_BASE}}/profile/conda/default"
        
        
        ${{ASSCOM2_BASE}}/asscom2 \
            --cores {wildcards.t} \
            --config \
                fofn={output.dir}/fofn.txt \
                prokka_rfam=false \
                output_directory={output.dir} \
            --until sequence_lengths assembly_stats prokka abricate mlst
        


    """
