
# Model

This directory contains the mathematical model and (inference framework)[INFERENCE.md].

### Parameters

#### Initial expansion and normal

\* Paremeter is calculated
Values don't represent those used.

| Parameter        | Value          | Description                   |
|------------------|----------------|-------------------------------|
| initSize         | 10000          | Initial size                  |
| birth_rate       | 0.05548397814  | Birth rate                    |
| death_rate       | 0.00548397814  | Death rate                    |
| mut_rate*         | 36.0           | Mutation rate                 |
| adv_mut_rate     | 0.0            | Advantageous mutation rate    |
| s_coef           | 0.0            | Selection coefficient         |
| max_t        | 10      | Maximum time                    |

#### Polyp Pop

| Parameter              | Value          | Description                        |
|------------------------|----------------|------------------------------------|
| num_seeds              | 1              | Number of seeds                    |
| final_pop_size         | 10000          | Final population size              |
| s_coef_polyp           | 0.0            | Selection coefficient for polyp    |
| adv_mut_rate_polyp     | 0.0            | Advantageous mutation rate for polyp |
| polyp_birth_rate       | 0.9            | Birth rate for polyp               |
| polyp_init_time        | 0.9            | Time of polyp initiation           |
| mut_rate_polyp         | 36.0           | Mutation rate for polyp            |
| max_t        | 10      | Maximum time (colectomy age)                    |

## IO and Options

| Parameter    | Value   | Description                     |
|--------------|---------|---------------------------------|
| outfile      | output  | Output path for file                     |
| seed         | 1234    | Random seed for simulations     |
| verbose      | 0       | Verbose output                  |
| mean_depth   | 100     | Mean sequencing depth           |
| sd_depth     | 10      | Standard deviation of sequencing depth |
| sample_size  | 10000   | Sample size of cells to take    |

### Output files and formats

The model outputs two files, `out.hdf5` and `muts.tsv`.

1. `out.hdf5` containing the model parameters, population sizes, mutID deliniation, and the the seed `genomes` used for the polyp initiation.

```
├── driver_muts
│   ├── fission (0)
│   ├── inutero (0)
│   └── polypfission (0)
├── mutId_delineator
│   └── delineator_muts (2)
├── params
│   ├── param_names (21)
│   └── values (21)
├── pops
│   ├── fission (21)
│   ├── inutero (10)
│   └── polypfission (21)
└── seed_genomes
    └── seed_genome_1 (411)
```

2. `muts.tsv` containing the sampled mutations that occurred in the population with the raw and 'sequenced' VAFs for each growth stage. Example:

```
mutation	raw_vaf	vaf	t	is_driver	seed_parent	group
5	0.11105	0.102803738317757	1	0	0	inutero
56	0.0605	0.07352941176470588	2	0	0	inutero
140	0.0164	0.017699115044247787	3	0	0	Normal
202	0.01075	0.011111111111111112	3	0	0	Normal
169	0.01265	0.010416666666666666	3	0	0	Normal
23921445	0.014511873350923483	0.0	13	0	1	Polyp
23921459	0.025065963060686015	0.03669724770642202	13	0	1	Polyp
23922455	0.0316622691292876	0.020202020202020204	15	0	1	Polyp
```

### Execution

The python wrapper uses a configuration file `conf.yaml` to run the model. The configuration file is used to set the parameters of the model, build output directory for model outputs and execute the model. An example of this is located in `TEST_Example.ipynb`.

```yaml
runtype: "full_model"  # Type of run

# Normal Pop
initSize: 10000  # Initial size
birth_rate: 0.05548397814  # Birth rate
death_rate: 0.00548397814  # Death rate
mut_rate: 12.0  # Mutation rate
adv_mut_rate: 0.0  # Advantageous mutation rate
s_coef: 0.0  # Selection coefficient

# Polyp Pop
num_seeds: 1  # Number of seeds
final_pop_size: 10000  # Final population size
s_coef_polyp: 0.0  # Selection coefficient for polyp
adv_mut_rate_polyp: 0.0  # Advantageous mutation rate for polyp
polyp_birth_rate: 0.9  # Birth rate for polyp
mut_rate_polyp: 12.0  # Mutation rate for polyp

# IO and Options
outfile: "output"  # Output file
seed: 1234  # Random seed for simulations
max_t: 10  # Maximum time
verbose: false  # Verbose output
mean_depth: 100  # Mean sequencing depth
sd_depth: 10  # Standard deviation of sequencing depth
sample_size: 10000  # Sample size of cells to take
```

#### Running the model without wrapper using `conf.yaml`:

```bash
$ julia ./model/run_model.jl --help
usage: run_model.jl [--runtype RUNTYPE] [--initSize INITSIZE]
                    [--birth_rate BIRTH_RATE]
                    [--death_rate DEATH_RATE] [--mut_rate MUT_RATE]
                    [--adv_mut_rate ADV_MUT_RATE] [--s_coef S_COEF]
                    [--num_seeds NUM_SEEDS]
                    [--final_pop_size FINAL_POP_SIZE]
                    [--s_coef_polyp S_COEF_POLYP]
                    [--adv_mut_rate_polyp ADV_MUT_RATE_POLYP]
                    [--polyp_birth_rate POLYP_BIRTH_RATE]
                    [--mut_rate_polyp MUT_RATE_POLYP]
                    [--outfile OUTFILE] [--seed SEED] [--max_t MAX_T]
                    [--verbose VERBOSE] [--mean_depth MEAN_DEPTH]
                    [--sd_depth SD_DEPTH] [--sample_size SAMPLE_SIZE]
                    [-h]

optional arguments:
  --runtype RUNTYPE     Type of run (default: "full_model")
  --initSize INITSIZE   Initial size (type: Int64, default: 10000)
  --birth_rate BIRTH_RATE
                        Birth rate (type: Float64, default: 0.055484)
  --death_rate DEATH_RATE
                        Death rate (type: Float64, default:
                        0.00548398)
  --mut_rate MUT_RATE   Mutation rate (type: Float64, default: 12.0)
  --adv_mut_rate ADV_MUT_RATE
                        Advantageous mutation rate (type: Float64,
                        default: 0.0)
  --s_coef S_COEF       Selection coefficient (type: Float64, default:
                        0.0)
  --num_seeds NUM_SEEDS
                        Number of seeds (type: Int64, default: 1)
  --final_pop_size FINAL_POP_SIZE
                        Final population size (type: Int64, default:
                        10000)
  --s_coef_polyp S_COEF_POLYP
                        Selection coefficient for polyp (type:
                        Float64, default: 0.0)
  --adv_mut_rate_polyp ADV_MUT_RATE_POLYP
                        Advantageous mutation rate for polyp (type:
                        Float64, default: 0.0)
  --polyp_birth_rate POLYP_BIRTH_RATE
                        Birth rate for polyp (type: Float64, default:
                        0.9)
  --mut_rate_polyp MUT_RATE_POLYP
                        Mutation rate for polyp (type: Float64,
                        default: 12.0)
  --outfile OUTFILE     Output file (default: "output")
  --seed SEED           Random seed for simulations (type: Int128,
                        default: 1234)
  --max_t MAX_T         Maximum time (type: Int64, default: 10)
  --verbose VERBOSE     Verbose output (type: Bool, default: false)
  --mean_depth MEAN_DEPTH
                        Mean sequencing depth (type: Int64, default:
                        100)
  --sd_depth SD_DEPTH   Standard deviation of sequencing depth (type:
                        Int64, default: 10)
  --sample_size SAMPLE_SIZE
                        Sample size of cells to take (type: Int64,
                        default: 10000)
  -h, --help            show this help message and exit
```

### Requirements

Developed using `Julia Version 1.9.4` and `Python 3.10.12`.

Julia packages:
```julia
using Pkg
Pkg.add(["FileIO", "CSV", "DataFrames", "HDF5", "Random", "Distributions", "ArgParse", "PoissonRandom", "StatsBase"])
```

Python packages:

Find these in `requirements.txt`


