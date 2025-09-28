# 2D-SPU
BRKGA-RLS algorithm for the Two-dimensional Rectangular Strip Packing Problem With Unloading Constraints Problem.

## Pre-requisites 🛠️

You must have C++17 and GNU g++ 11.4.0

## How to run ⚙️

```console
$ ./brkgarls -f instance_path -s seed + [BRKGA-RLS-PARAMS]
```

or

```console
$ ./run.sh
```

Available flags:

* **-f <file>**: Specifies the path of the input instance. 📂 [required]
* **-s <string>**: Sets the seed for generating random numbers. 🌱 [required]
* **-t <number>**: Sets the time limit for execution (defaults to 60s). ⏳
* **-d**: Activates debug mode, the output is a file you can use as input in the `debug.py` script. 🐞
* **-h**: Shows this help menu. 📖

BRKGA-RLS-PARAMS [required]:
* **-p <number>**: Number of individuals in each population. 👥
* **-g <number>**: Frequency (in generations) of elite migration between populations. 🔄
* **-n <number>**: Number of elite individuals transferred between populations. 🚚
* **-e <float>**: Fraction of the population formed by elite. ⭐
* **-m <float>**: Fraction of the population formed by mutants. 🧬
* **-o <float>**: Probability of inheriting an allele from the elite parent (*rho_e*). 🎲
* **-k <number>**: Number of independent populations (*K*). 🌍
* **-x <float>**: Fraction of the initial population from naive solutions (*gamma_1*). 🔹
* **-y <float>**: Fraction of the initial population from naive solutions (*gamma_2*). 🔸
* **-z <float>**: Fraction of the initial population from naive solutions (*gamma_3*). ▪️
* **-i <number>**: Swap window size used in initial solution generation (*w_i*). 🪟
* **-a <number>**: Number of RLS iterations (*i_RLS*). 🔁
* **-l <number>**: Swap window size used in RLS (*w_RLS*). 🪟
* **-v <float>**: Probability of applying RLS to an elite individual (*rho_RLS*). 🎯

## How to debug 🐞

An auxiliary Python tool was developed to help debug the strip packing layout solution. It plots the arranged items on the strip and their maximum width to visualize the generated solutions.

**Visualization Details:**
* Numbers on the items indicate their position in the item vector.
* Items from the same client share the same color.
* The color palette maps client classes: **Higher classes** are closer to **red**, and **lower classes** are closer to **violet**.

```console
$ ./debug.py solution_file_path
```

## Metrics

```console
$ ./metrics.py package_name_in_logs_folder
```
