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

dengue2 = myDir + '/dengue2'
//myDir.eachFileMatch ~/.*_1.fastq.gz/, {L001R1Lst << it.name}
dengue2.eachFileMatch ~/.*_1.fastq.gz/, {L001R1Lst << it.name}

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

include { quality2 } from './modules/quality2.nf'
include { nofrag2 } from './modules/nofrag2.nf'
include { frag2 } from './modules/frag2.nf'
include { primer2 } from './modules/primer2.nf'
include { assembly2 } from './modules/assembly2.nf'
include { pystats2 } from './modules/pystats2.nf'

workflow {
    if("${params.frag}" == "frag"){
       quality2(A) | frag2 | primer2 | assembly2 | pystats2 | view
    }
    else{
       quality2(A) | nofrag2 | primer2 | assembly2 | pystats2 | view
    }
}


