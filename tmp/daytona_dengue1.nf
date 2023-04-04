#!/usr/bin/env nextflow

/*
Note:
Before running the script, please set the parameters in the config file params.yaml
*/

//Step1:input data files
nextflow.enable.dsl=2
def L001R1Lst = []
def sampleNames = []
myDir = file("$params.input")

dengue1 = myDir + '/dengue1'
//myDir.eachFileMatch ~/.*_1.fastq.gz/, {L001R1Lst << it.name}
dengue1.eachFileMatch ~/.*_1.fastq.gz/, {L001R1Lst << it.name}

L001R1Lst.sort()
L001R1Lst.each{
   def x = it.minus("_1.fastq.gz")
     //println x
   sampleNames.add(x)
}
//println L001R1Lst
//println sampleNames


//Step2: process the inputed data
A = Channel.fromList(sampleNames)
//A.view()

include { quality1 } from './modules/quality1.nf'
include { nofrag1 } from './modules/nofrag1.nf'
include { frag1 } from './modules/frag1.nf'
include { primer1 } from './modules/primer1.nf'
include { assembly1 } from './modules/assembly1.nf'
include { pystats1 } from './modules/pystats1.nf'

workflow {
    if("${params.frag}" == "frag"){
       quality1(A) | frag1 | primer1 | assembly1 | pystats1 | view
    }
    else{
       quality(A) | nofrag1 | primer1 | assembly1 | pystats1 | view
    }
}


