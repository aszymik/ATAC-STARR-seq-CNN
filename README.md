## Using convolution neural networks (CNN) for better understanding of human genome

This project constitutes a convolutional neural network for predicting regulatory regions in a multi-class classification task.
It is an extension of https://github.com/marnifora/magisterka repository.

### Versions
Python 3.10

All packages versions can be found in *requirements.txt* file.


### Example models

Example of trained networks are in */data/custom40* and */data/alt1* directories.


### Example data

Sample data containing longer sequences from ChIP-seq experiments are located in *data/test_fasta* folder, while data with shorter transcriptionally active regions from ATAC-STARR-seq are located in *data/test_shorter_seq* directory.


### Training an existing network on new samples

To re-train a network on newly available data, you can use the *train_new.py* script.

The simplest way to invoke is as follows:

```
 python3 train_new.py 
```

It will then run an analysis based on the data in the *data/test_fasta* and the model from *data/alt1* directory.

If you this runs properly you can start your own analysis by changing some parameters. Let us start by looking at the options:

```
$ python3 train_new.py --help
usage: train_new.py [-h] [--model NAME] [-p DIR] [-x PREFIX] [-o DIR] [--seed NUMBER] [--optimizer NAME] [--loss_fn NAME] [--batch_size INT] [--num_workers INT] [--num_epochs INT]
                    [--acc_threshold FLOAT] [--learning_rate FLOAT] [--no_adjust_lr] [--seq_len INT] [--dropout-conv FLOAT] [--dropout-fc FLOAT] [--weight-decay FLOAT] [--momentum FLOAT]

Train network based on given data

options:
 -h, --help            show this help message and exit
  --model NAME          File with the model weights to load before training
  --namespace NAME      The namespace for this run of training
  -p DIR, --path DIR    Working directory.
  -x PREFIX, --prefix PREFIX
                        file_prefix.
  -o DIR, --output DIR  Output directory, default: test_output
  --seed NUMBER         Set random seed, default: 0
  --optimizer NAME      optimization algorithm to use for training the network, default = RMSprop
  --loss_fn NAME        loss function for training the network, default = CrossEntropyLoss
  --batch_size INT      size of the batch, default: 64
  --num_workers INT     how many subprocesses to use for data loading, default: 4
  --num_epochs INT      maximum number of epochs to run, default: 300
  --acc_threshold FLOAT
                        threshold of the validation accuracy - if gained training process stops, default: 0.9
  --learning_rate FLOAT
                        initial learning rate, default: 0.01
  --no_adjust_lr        no reduction of learning rate during training, default: False
  --seq_len INT         Length of the input sequences to the network, default: 2000
  --dropout-conv FLOAT  Dropout of convolutional layers, default value is 0.2
  --dropout-fc FLOAT    Dropout of fully-connected layers, default value is 0.5
  --weight-decay FLOAT  Weight decay, default value is 0.0001
  --momentum FLOAT      Momentum, default value is 0.1
  --no_random_noise     Not replacing Ns in the sequence with random nucleotides, default: False
```

If we just want to use a different fasta files, without modifying the training params, you can simply run:

```
python3 --path data/my_fasta --prefix fasta_file -output data/my_output --namespace MY-RUN-1
```

for this we assume that you have 4 fasta files in the *data/my_fasta* folder and that the *data/my_output* folder is created. 



### Calculating the outputs for a full fasta file

If you just want to obtain the output of a network processing a fasta file, you can use the following command:

```
python3 process_fasta.py -m data/alt1/alt1_last.model \\
        -i data/test_fasta/test_pa.fa -o data/test_output/test_pa_out
```

You need to specify the network model, input fasta file and output text file name. It can take a while, but you will end up with a file with the names of all seuqnces and the values of all ouput neurons on them.



### Calculating integrated gradients

To calculate integrads based on example model and set of sequences just run:

```
python3 calculate_integrads.py \
        --model custom40_last.model \
        --seq extreme_custom40_train_1.fasta \
        --baseline CHOSEN-BASELINE
```
CHOSEN-BASELINE depends on what baseline you want to use for calculating 
integrated gradients (see: https://arxiv.org/abs/1703.01365 for details 
of the method), select one of the options:
- *zeros* - use zeros array as baseline
- *fixed* - use the same set of random sequences as the baseline for each 
sequence
- *random* - use different random set of sequences as the baseline for each 
sequence
- *test-balanced-8-3_baseline.npy* - use pre-calculated balanced baseline 
(each 3 nucleotides in the given position occur exactly once across all 64. baseline sequences)

As the output new directory called *integrads_NETWORK_SEQUENCES_BASELINE_TRIALS-STEPS* is created
(*integrads_custom40_extreme-custom40-train-1_CHOSEN-BASELINE_10-50* if the default data was used). 
Inside there are result files:
- integrads_all.npy - numpy array with calculated gradients
- params.txt - file with parameters of the analysis

### Plotting seqlogos

To plot seqlogos based on the calculated integrads run:
```
python3 plot_seqlogo.py \
integrads_custom40_extreme-custom40-train-1_CHOSEN-BASELINE_10-50/integrads_all.npy \
--global_ylim \
--one \
--clip NUM
```
Options *global_ylim*, *one* and *clip* are optional:
- *global_ylim* - 
set the same y axis range for all sequences from the given array
- *one* - plot only one nucleotide in one position 
(the one from the original sequence)
- *clip* - subset of nucleotides to plot: +-NUM from the middle 
of the sequences, by default NUM=100

As the output new directory with plots is created.


## Analysing shorter sequences

To prepare network input from example shorter sequences, run:

```
python3 pad_fasta_with_Ns.py  data/test_shorter_seq/shorter_seq_pa.fa  OUTPUT_PATH
```
You can specify the minimum sequence length to be included as well as the target sequence length using the `--min_length` and `--output_length` arguments.


To truncate longer sequences to a given length and pad them with Ns, run:
```
python3 create_shorter_seq.py  data/test_shorter_seq/2000bp_pi.fa  OUTPUT_PATH \
        -l DESIRED_LENGTH
```

To shorten longer sequences heterogeneously, based on the distribution of other sequences set, run the `change_seq_length_distributions.py` script.

Its arguments include:
```
  --reference           FASTA file paths of which distribution is to be mimicked (at least one is required)
  --input INPUT         Input FASTA file path (an example is *data/test_shorter_seq/2000bp_pi.fa*)
  --output OUTPUT       Output FASTA file path
  --length LENGTH       Target output length for sequences
  --min_length          Minimum length for sequences to be included
  --subset_size         Size of the subset to create
  --exact               If true, samples from exact lengths from distribution, otherwise from calculated
                        inverse cumulative distribution function; default: True
```

Example run:
```
python3 change_seq_length_distributions.py \
        -r data/test_shorter_seq/shorter_seq_pa.fa data/test_shorter_seq/shorter_seq_na.fa \
        -i data/test_shorter_seq/2000bp_pi.fa \
        -o data/my_output.fa
```

