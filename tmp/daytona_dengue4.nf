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

dengue4 = myDir + '/dengue4'
//myDir.eachFileMatch ~/.*_1.fastq.gz/, {L001R1Lst << it.name}
dengue4.eachFileMatch ~/.*_1.fastq.gz/, {L001R1Lst << it.name}

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

include { quality4 } from './modules/quality4.nf'
include { nofrag4 } from './modules/nofrag4.nf'
include { frag4 } from './modules/frag4.nf'
include { primer4 } from './modules/primer4.nf'
include { assembly4 } from './modules/assembly4.nf'
include { pystats4 } from './modules/pystats4.nf'

workflow {
    if("${params.frag}" == "frag"){
       quality4(A) | frag4 | primer4 | assembly4 | pystats4 | view
    }
    else{
       quality4(A) | nofrag4 | primer4 | assembly4 | pystats4 | view
    }
}


