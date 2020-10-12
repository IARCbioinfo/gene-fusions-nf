#!/usr/bin/env nextflow

params.in = "$baseDir/data/sample.fa"

/*
 * split a fasta file in multiple files
 */
process splitSequences {

    input:
    path 'input.fa' from params.in

    output:
    path 'seq_*' into records

    """
    awk '/^>/{f="seq_"++d} {print > f}' < input.fa
    arriba --version
    star --version
    """
}

/*
 * Simple reverse the sequences
 */
process reverse {

    input:
    path x from records

    output:
    stdout into result

    """
    cat $x | rev
    """
}

/*
 * print the channel content
 */
result.subscribe { println it }
