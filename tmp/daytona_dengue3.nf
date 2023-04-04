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

dengue3 = myDir + '/dengue3'
//myDir.eachFileMatch ~/.*_1.fastq.gz/, {L001R1Lst << it.name}
dengue3.eachFileMatch ~/.*_1.fastq.gz/, {L001R1Lst << it.name}

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

include { quality3 } from './modules/quality3.nf'
include { nofrag3 } from './modules/nofrag3.nf'
include { frag3 } from './modules/frag3.nf'
include { primer3 } from './modules/primer3.nf'
include { assembly3 } from './modules/assembly3.nf'
include { pystats3 } from './modules/pystats3.nf'

workflow {
    if("${params.frag}" == "frag"){
       quality3(A) | frag3 | primer3 | assembly3 | pystats3 | view
    }
    else{
       quality3(A) | nofrag3 | primer3 | assembly3 | pystats3 | view
    }
}


