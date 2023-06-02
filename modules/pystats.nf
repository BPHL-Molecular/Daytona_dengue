process pystats {
    input:
        val samplename
    output:
        stdout
        //val mypath
        //path "pyoutputs.txt", emit: pyoutputs
        
    $/
    #!/usr/bin/env python3
    import subprocess
    import re
    
    if re.match("SER1_", "${samplename}"):
       print("SER1, ${samplename}")
       
       filepath1 = "${params.output}"+"/dengue1/"+"${samplename}"+"/alignment/"+"${samplename}"+".coverage.txt"
       #print(filepath1)
       with open(filepath1, 'r') as cov_report:
           header = cov_report.readline()
           header = header.rstrip()
           stats = cov_report.readline()
           stats = stats.rstrip()
           stats = stats.split()
           ref_name = stats[0]
           #print(ref_name)
           start = stats[1]
           end = stats[2]
           reads_mapped = stats[3]
           cov_bases = stats[4]
           cov = stats[5]
           depth = stats[6]
           baseq = stats[7]
           #print(reads_mapped)
           mapq = stats[8]
        
       #Get number of raw reads
       proc_1 = subprocess.run('zcat ' + "${params.output}"+"/dengue1/"+"${samplename}"+"/"+"${samplename}" + '_1_humanclean.fastq.gz | wc -l', shell=True, capture_output=True, text=True, check=True)
       wc_out_1 = proc_1.stdout.rstrip()
       reads_1 = int(wc_out_1) / 4
       proc_2 = subprocess.run('zcat ' + "${params.output}"+"/dengue1/"+"${samplename}"+"/"+"${samplename}" + '_2_humanclean.fastq.gz | wc -l', shell=True, capture_output=True, text=True, check=True)
       wc_out_2 = proc_2.stdout.rstrip()
       reads_2 = int(wc_out_2) / 4
       raw_reads = reads_1 + reads_2
       raw_reads = int(raw_reads)

       #Get number of clean reads
       proc_c1x = subprocess.run('zcat ' + "${params.output}"+"/dengue1/"+"${samplename}"+"/"+"${samplename}" + '_1.fq.gz | wc -l', shell=True, capture_output=True, text=True, check=True)
       wc_out_c1x = proc_c1x.stdout.rstrip()
       reads_c1x = int(wc_out_c1x) / 4
       proc_c2x = subprocess.run('zcat ' + "${params.output}"+"/dengue1/"+"${samplename}"+"/"+"${samplename}" + '_2.fq.gz | wc -l', shell=True, capture_output=True, text=True, check=True)
       wc_out_c2x = proc_c2x.stdout.rstrip()
       reads_c2x = int(wc_out_c2x) / 4
       clean_reads = reads_c1x + reads_c2x
       clean_reads = int(clean_reads)
       #print(clean_reads)
    
       #Get percentage of mapped reads/clean reads
       percent_map = "%0.4f"%((int(reads_mapped)/int(clean_reads))*100)
       #print(percent_map)
    
       #Gather QC metrics for consensus assembly
       filepath2 = "${params.output}"+"/dengue1/"+"${samplename}"+"/assembly/"+"${samplename}"+".consensus.fa"
       with open(filepath2, 'r') as assem:
           header = assem.readline()
           header = header.rstrip()
           bases = assem.readline()
           bases = bases.rstrip()
           num_bases = len(bases)
           ns = bases.count('N')
           called = num_bases - ns
           pg = "%0.4f"%((called/int(end))*100)
           #print(called)
           #print(end)
           #print(pg)
       #Rename header in fasta to just sample name
       subprocess.run("sed -i \'s/^>.*/>"+"${samplename}"+"/\' "+filepath2, shell=True, check=True)
    
       #QC flag
       pg_flag = ''
       dp_flag = ''
       qc_flag = ''
       if float(pg) < 79.5:
           pg_flag = 'FAIL: Percent genome < 80%'
           qc_flag = qc_flag + pg_flag
       else:
           if float(depth) < 100:
              dp_flag = 'FAIL: Mean read depth < 100x'
              qc_flag = qc_flag + dp_flag
           if qc_flag == '':
              qc_flag = qc_flag + 'PASS'
       #print(qc_flag)
    
       idname = "${samplename}".lstrip("SER1_")

       if qc_flag == 'PASS':
             
           subprocess.run("cp "+filepath2+" "+"${params.output}"+"/dengue1/assemblies/"+idname+".consensus.fa", shell=True, check=True)   
           subprocess.run('cp ' + "${params.output}"+"/dengue1/"+"${samplename}" + '/variants/'+"${samplename}" + '.variants.tsv ' + "${params.output}"+'/dengue1/variants/'+idname+'.variants.tsv', shell=True, check=True)
    
           
    
           #Run VADR
           out_log = open("${params.output}"+"/dengue1/"+"${samplename}"+"/"+"${samplename}"+'.out', 'w')
           err_log = open("${params.output}"+"/dengue1/"+"${samplename}"+"/"+"${samplename}"+'.err', 'w')
           subprocess.run("singularity exec -B "+"${params.output}"+"/dengue1/"+"${samplename}"+"/assembly"+":/data docker://staphb/vadr:1.5.1 /opt/vadr/vadr/miniscripts/fasta-trim-terminal-ambigs.pl --minlen 50 --maxlen 30000 /data/" + "${samplename}"+".consensus.fa > " + "${params.output}"+"/dengue1/"+"${samplename}"+"/assembly/"+"${samplename}"+".trimmed.fasta", shell=True, stdout=out_log, stderr=err_log, check=True)

           subprocess.run("singularity exec -B "+"${params.output}"+"/dengue1/"+"${samplename}"+"/assembly:/data docker://staphb/vadr:1.5.1 v-annotate.pl --split --cpu 8 --group Dengue --nomisc --noprotid --mkey flavi --mdir /opt/vadr/vadr-models-flavi/ --noseqnamemax /data/"+"${samplename}"+".trimmed.fasta /data/"+"vadr_results", shell=True, stdout=out_log, stderr=err_log, check=True)


           #Parse through VADR outputs to get PASS or REVIEW flag
           vadr_flag = ''
           with open("${params.output}"+"/dengue1/"+"${samplename}"+"/assembly/vadr_results/vadr_results.vadr.pass.list", 'r') as p_list:
               result = p_list.readline()
               result = result.rstrip()
               if result == "${samplename}":
                   vadr_flag = 'PASS'

           with open("${params.output}"+"/dengue1/"+"${samplename}"+"/assembly/vadr_results/vadr_results.vadr.fail.list", 'r') as f_list:
               result = f_list.readline()
               result = result.rstrip()
               if result == "${samplename}":
                   vadr_flag = 'REVIEW'

           #Copy VADR error report to main analysis folder for easier review
           if vadr_flag == 'REVIEW':
               subprocess.run("cp " + "${params.output}"+"/dengue1/"+"${samplename}"+"/assembly/vadr_results/vadr_results.vadr.alt.list "+"${params.output}"+"/dengue1/vadr_error_reports/"+"${samplename}"+".vadr.alt.list", shell=True, check=True)

           out_log.close()
           err_log.close()
           subprocess.run("mv " + "${params.output}"+"/dengue1/"+"${samplename}"+"/assembly/vadr_results "+"${params.output}"+"/dengue1/"+"${samplename}"+"/assembly/"+"${samplename}"+"_vadr_results", shell=True, check=True)

       else:
           vadr_flag = 'NA'
           lineage = 'NA'
           sotc_out = 'NA'
        
       with open("${params.output}"+"/dengue1/"+"${samplename}"+"/report.txt", 'w') as report:
           header = ['sampleID', 'reference', 'start', 'end', 'num_raw_reads', 'num_clean_reads', 'num_mapped_reads', 'percent_mapped_clean_reads', 'cov_bases_mapped', 'percent_genome_cov_map', 'mean_depth', 'mean_base_qual', 'mean_map_qual', 'assembly_length', 'numN', 'percent_ref_genome_cov', 'VADR_flag', 'QC_flag']
           report.write('\t'.join(map(str,header)) + '\n')
           results = [idname, ref_name, start, end, raw_reads, clean_reads, reads_mapped, percent_map, cov_bases, cov, depth, baseq, mapq, num_bases, ns, pg, vadr_flag, qc_flag]
           report.write('\t'.join(map(str,results)) + '\n')
       

    elif re.match("SER2_", "${samplename}"):
       print("SER2, ${samplename}")
       
       filepath1 = "${params.output}"+"/dengue2/"+"${samplename}"+"/alignment/"+"${samplename}"+".coverage.txt"
       #print(filepath1)
       with open(filepath1, 'r') as cov_report:
           header = cov_report.readline()
           header = header.rstrip()
           stats = cov_report.readline()
           stats = stats.rstrip()
           stats = stats.split()
           ref_name = stats[0]
           #print(ref_name)
           start = stats[1]
           end = stats[2]
           reads_mapped = stats[3]
           cov_bases = stats[4]
           cov = stats[5]
           depth = stats[6]
           baseq = stats[7]
           #print(reads_mapped)
           mapq = stats[8]
        
       #Get number of raw reads
       proc_1 = subprocess.run('zcat ' + "${params.output}"+"/dengue2/"+"${samplename}"+"/"+"${samplename}" + '_1_humanclean.fastq.gz | wc -l', shell=True, capture_output=True, text=True, check=True)
       wc_out_1 = proc_1.stdout.rstrip()
       reads_1 = int(wc_out_1) / 4
       proc_2 = subprocess.run('zcat ' + "${params.output}"+"/dengue2/"+"${samplename}"+"/"+"${samplename}" + '_2_humanclean.fastq.gz | wc -l', shell=True, capture_output=True, text=True, check=True)
       wc_out_2 = proc_2.stdout.rstrip()
       reads_2 = int(wc_out_2) / 4
       raw_reads = reads_1 + reads_2
       raw_reads = int(raw_reads)

       #Get number of clean reads
       proc_c1x = subprocess.run('zcat ' + "${params.output}"+"/dengue2/"+"${samplename}"+"/"+"${samplename}" + '_1.fq.gz | wc -l', shell=True, capture_output=True, text=True, check=True)
       wc_out_c1x = proc_c1x.stdout.rstrip()
       reads_c1x = int(wc_out_c1x) / 4
       proc_c2x = subprocess.run('zcat ' + "${params.output}"+"/dengue2/"+"${samplename}"+"/"+"${samplename}" + '_2.fq.gz | wc -l', shell=True, capture_output=True, text=True, check=True)
       wc_out_c2x = proc_c2x.stdout.rstrip()
       reads_c2x = int(wc_out_c2x) / 4
       clean_reads = reads_c1x + reads_c2x
       clean_reads = int(clean_reads)
       #print(clean_reads)
    
       #Get percentage of mapped reads/clean reads
       percent_map = "%0.4f"%((int(reads_mapped)/int(clean_reads))*100)
       #print(percent_map)
    
       #Gather QC metrics for consensus assembly
       filepath2 = "${params.output}"+"/dengue2/"+"${samplename}"+"/assembly/"+"${samplename}"+".consensus.fa"
       with open(filepath2, 'r') as assem:
           header = assem.readline()
           header = header.rstrip()
           bases = assem.readline()
           bases = bases.rstrip()
           num_bases = len(bases)
           ns = bases.count('N')
           called = num_bases - ns
           pg = "%0.4f"%((called/int(end))*100)
           #print(called)
           #print(end)
           #print(pg)
       #Rename header in fasta to just sample name
       subprocess.run("sed -i \'s/^>.*/>"+"${samplename}"+"/\' "+filepath2, shell=True, check=True)
    
       #QC flag
       pg_flag = ''
       dp_flag = ''
       qc_flag = ''
       if float(pg) < 79.5:
           pg_flag = 'FAIL: Percent genome < 80%'
           qc_flag = qc_flag + pg_flag
       else:
           if float(depth) < 100:
              dp_flag = 'FAIL: Mean read depth < 100x'
              qc_flag = qc_flag + dp_flag
           if qc_flag == '':
              qc_flag = qc_flag + 'PASS'
       #print(qc_flag)
       
       idname = "${samplename}".lstrip("SER2_")
       if qc_flag == 'PASS':          
           subprocess.run("cp "+filepath2+" "+"${params.output}"+"/dengue2/assemblies/"+idname+".consensus.fa", shell=True, check=True)   
           subprocess.run('cp ' + "${params.output}"+"/dengue2/"+"${samplename}" + '/variants/'+"${samplename}" + '.variants.tsv ' + "${params.output}"+'/dengue2/variants/'+idname+'.variants.tsv', shell=True, check=True)
    
           #Run VADR
           out_log = open("${params.output}"+"/dengue2/"+"${samplename}"+"/"+"${samplename}"+'.out', 'w')
           err_log = open("${params.output}"+"/dengue2/"+"${samplename}"+"/"+"${samplename}"+'.err', 'w')
           subprocess.run("singularity exec -B "+"${params.output}"+"/dengue2/"+"${samplename}"+"/assembly"+":/data docker://staphb/vadr:1.5.1 /opt/vadr/vadr/miniscripts/fasta-trim-terminal-ambigs.pl --minlen 50 --maxlen 30000 /data/" + "${samplename}"+".consensus.fa > " + "${params.output}"+"/dengue2/"+"${samplename}"+"/assembly/"+"${samplename}"+".trimmed.fasta", shell=True, stdout=out_log, stderr=err_log, check=True)

           subprocess.run("singularity exec -B "+"${params.output}"+"/dengue2/"+"${samplename}"+"/assembly:/data docker://staphb/vadr:1.5.1 v-annotate.pl --split --cpu 8 --group Dengue --nomisc --noprotid --mkey flavi --mdir /opt/vadr/vadr-models-flavi/ --noseqnamemax /data/"+"${samplename}"+".trimmed.fasta /data/"+"vadr_results", shell=True, stdout=out_log, stderr=err_log, check=True)


           #Parse through VADR outputs to get PASS or REVIEW flag
           vadr_flag = ''
           with open("${params.output}"+"/dengue2/"+"${samplename}"+"/assembly/vadr_results/vadr_results.vadr.pass.list", 'r') as p_list:
               result = p_list.readline()
               result = result.rstrip()
               if result == "${samplename}":
                   vadr_flag = 'PASS'

           with open("${params.output}"+"/dengue2/"+"${samplename}"+"/assembly/vadr_results/vadr_results.vadr.fail.list", 'r') as f_list:
               result = f_list.readline()
               result = result.rstrip()
               if result == "${samplename}":
                   vadr_flag = 'REVIEW'

           #Copy VADR error report to main analysis folder for easier review
           if vadr_flag == 'REVIEW':
               subprocess.run("cp " + "${params.output}"+"/dengue2/"+"${samplename}"+"/assembly/vadr_results/vadr_results.vadr.alt.list "+"${params.output}"+"/dengue2/vadr_error_reports/"+"${samplename}"+".vadr.alt.list", shell=True, check=True)

           out_log.close()
           err_log.close()
           subprocess.run("mv " + "${params.output}"+"/dengue2/"+"${samplename}"+"/assembly/vadr_results "+"${params.output}"+"/dengue2/"+"${samplename}"+"/assembly/"+"${samplename}"+"_vadr_results", shell=True, check=True)

       else:
           vadr_flag = 'NA'
           lineage = 'NA'
           sotc_out = 'NA'
        
       with open("${params.output}"+"/dengue2/"+"${samplename}"+"/report.txt", 'w') as report:
           header = ['sampleID', 'reference', 'start', 'end', 'num_raw_reads', 'num_clean_reads', 'num_mapped_reads', 'percent_mapped_clean_reads', 'cov_bases_mapped', 'percent_genome_cov_map', 'mean_depth', 'mean_base_qual', 'mean_map_qual', 'assembly_length', 'numN', 'percent_ref_genome_cov', 'VADR_flag', 'QC_flag']
           report.write('\t'.join(map(str,header)) + '\n')
           results = [idname, ref_name, start, end, raw_reads, clean_reads, reads_mapped, percent_map, cov_bases, cov, depth, baseq, mapq, num_bases, ns, pg, vadr_flag, qc_flag]
           report.write('\t'.join(map(str,results)) + '\n')


    elif re.match("SER3_", "${samplename}"):
       print("SER3, ${samplename}")
       
       filepath1 = "${params.output}"+"/dengue3/"+"${samplename}"+"/alignment/"+"${samplename}"+".coverage.txt"
       #print(filepath1)
       with open(filepath1, 'r') as cov_report:
           header = cov_report.readline()
           header = header.rstrip()
           stats = cov_report.readline()
           stats = stats.rstrip()
           stats = stats.split()
           ref_name = stats[0]
           #print(ref_name)
           start = stats[1]
           end = stats[2]
           reads_mapped = stats[3]
           cov_bases = stats[4]
           cov = stats[5]
           depth = stats[6]
           baseq = stats[7]
           #print(reads_mapped)
           mapq = stats[8]
        
       #Get number of raw reads
       proc_1 = subprocess.run('zcat ' + "${params.output}"+"/dengue3/"+"${samplename}"+"/"+"${samplename}" + '_1_humanclean.fastq.gz | wc -l', shell=True, capture_output=True, text=True, check=True)
       wc_out_1 = proc_1.stdout.rstrip()
       reads_1 = int(wc_out_1) / 4
       proc_2 = subprocess.run('zcat ' + "${params.output}"+"/dengue3/"+"${samplename}"+"/"+"${samplename}" + '_2_humanclean.fastq.gz | wc -l', shell=True, capture_output=True, text=True, check=True)
       wc_out_2 = proc_2.stdout.rstrip()
       reads_2 = int(wc_out_2) / 4
       raw_reads = reads_1 + reads_2
       raw_reads = int(raw_reads)

       #Get number of clean reads
       proc_c1x = subprocess.run('zcat ' + "${params.output}"+"/dengue3/"+"${samplename}"+"/"+"${samplename}" + '_1.fq.gz | wc -l', shell=True, capture_output=True, text=True, check=True)
       wc_out_c1x = proc_c1x.stdout.rstrip()
       reads_c1x = int(wc_out_c1x) / 4
       proc_c2x = subprocess.run('zcat ' + "${params.output}"+"/dengue3/"+"${samplename}"+"/"+"${samplename}" + '_2.fq.gz | wc -l', shell=True, capture_output=True, text=True, check=True)
       wc_out_c2x = proc_c2x.stdout.rstrip()
       reads_c2x = int(wc_out_c2x) / 4
       clean_reads = reads_c1x + reads_c2x
       clean_reads = int(clean_reads)
       #print(clean_reads)
    
       #Get percentage of mapped reads/clean reads
       percent_map = "%0.4f"%((int(reads_mapped)/int(clean_reads))*100)
       #print(percent_map)
    
       #Gather QC metrics for consensus assembly
       filepath2 = "${params.output}"+"/dengue3/"+"${samplename}"+"/assembly/"+"${samplename}"+".consensus.fa"
       with open(filepath2, 'r') as assem:
           header = assem.readline()
           header = header.rstrip()
           bases = assem.readline()
           bases = bases.rstrip()
           num_bases = len(bases)
           ns = bases.count('N')
           called = num_bases - ns
           pg = "%0.4f"%((called/int(end))*100)
           #print(called)
           #print(end)
           #print(pg)
       #Rename header in fasta to just sample name
       subprocess.run("sed -i \'s/^>.*/>"+"${samplename}"+"/\' "+filepath2, shell=True, check=True)
    
       #QC flag
       pg_flag = ''
       dp_flag = ''
       qc_flag = ''
       if float(pg) < 79.5:
           pg_flag = 'FAIL: Percent genome < 80%'
           qc_flag = qc_flag + pg_flag
       else:
           if float(depth) < 100:
              dp_flag = 'FAIL: Mean read depth < 100x'
              qc_flag = qc_flag + dp_flag
           if qc_flag == '':
              qc_flag = qc_flag + 'PASS'
       #print(qc_flag)
       
       idname = "${samplename}".lstrip("SER3_")
       if qc_flag == 'PASS':           
           subprocess.run("cp "+filepath2+" "+"${params.output}"+"/dengue3/assemblies/"+idname+".consensus.fa", shell=True, check=True)   
           subprocess.run('cp ' + "${params.output}"+"/dengue3/"+"${samplename}" + '/variants/'+"${samplename}" + '.variants.tsv ' + "${params.output}"+'/dengue3/variants/'+idname+'.variants.tsv', shell=True, check=True)
    
           #Run VADR
           out_log = open("${params.output}"+"/dengue3/"+"${samplename}"+"/"+"${samplename}"+'.out', 'w')
           err_log = open("${params.output}"+"/dengue3/"+"${samplename}"+"/"+"${samplename}"+'.err', 'w')
           subprocess.run("singularity exec -B "+"${params.output}"+"/dengue3/"+"${samplename}"+"/assembly"+":/data docker://staphb/vadr:1.5.1 /opt/vadr/vadr/miniscripts/fasta-trim-terminal-ambigs.pl --minlen 50 --maxlen 30000 /data/" + "${samplename}"+".consensus.fa > " + "${params.output}"+"/dengue3/"+"${samplename}"+"/assembly/"+"${samplename}"+".trimmed.fasta", shell=True, stdout=out_log, stderr=err_log, check=True)

           subprocess.run("singularity exec -B "+"${params.output}"+"/dengue3/"+"${samplename}"+"/assembly:/data docker://staphb/vadr:1.5.1 v-annotate.pl --split --cpu 8 --group Dengue --nomisc --noprotid --mkey flavi --mdir /opt/vadr/vadr-models-flavi/ --noseqnamemax /data/"+"${samplename}"+".trimmed.fasta /data/"+"vadr_results", shell=True, stdout=out_log, stderr=err_log, check=True)


           #Parse through VADR outputs to get PASS or REVIEW flag
           vadr_flag = ''
           with open("${params.output}"+"/dengue3/"+"${samplename}"+"/assembly/vadr_results/vadr_results.vadr.pass.list", 'r') as p_list:
               result = p_list.readline()
               result = result.rstrip()
               if result == "${samplename}":
                   vadr_flag = 'PASS'

           with open("${params.output}"+"/dengue3/"+"${samplename}"+"/assembly/vadr_results/vadr_results.vadr.fail.list", 'r') as f_list:
               result = f_list.readline()
               result = result.rstrip()
               if result == "${samplename}":
                   vadr_flag = 'REVIEW'

           #Copy VADR error report to main analysis folder for easier review
           if vadr_flag == 'REVIEW':
               subprocess.run("cp " + "${params.output}"+"/dengue3/"+"${samplename}"+"/assembly/vadr_results/vadr_results.vadr.alt.list "+"${params.output}"+"/dengue3/vadr_error_reports/"+"${samplename}"+".vadr.alt.list", shell=True, check=True)

           out_log.close()
           err_log.close()
           subprocess.run("mv " + "${params.output}"+"/dengue3/"+"${samplename}"+"/assembly/vadr_results "+"${params.output}"+"/dengue3/"+"${samplename}"+"/assembly/"+"${samplename}"+"_vadr_results", shell=True, check=True)

       else:
           vadr_flag = 'NA'
           lineage = 'NA'
           sotc_out = 'NA'
        
       with open("${params.output}"+"/dengue3/"+"${samplename}"+"/report.txt", 'w') as report:
           header = ['sampleID', 'reference', 'start', 'end', 'num_raw_reads', 'num_clean_reads', 'num_mapped_reads', 'percent_mapped_clean_reads', 'cov_bases_mapped', 'percent_genome_cov_map', 'mean_depth', 'mean_base_qual', 'mean_map_qual', 'assembly_length', 'numN', 'percent_ref_genome_cov', 'VADR_flag', 'QC_flag']
           report.write('\t'.join(map(str,header)) + '\n')
           results = [idname, ref_name, start, end, raw_reads, clean_reads, reads_mapped, percent_map, cov_bases, cov, depth, baseq, mapq, num_bases, ns, pg, vadr_flag, qc_flag]
           report.write('\t'.join(map(str,results)) + '\n')


    elif re.match("SER4_", "${samplename}"):
       print("SER4, ${samplename}")
       
       filepath1 = "${params.output}"+"/dengue4/"+"${samplename}"+"/alignment/"+"${samplename}"+".coverage.txt"
       #print(filepath1)
       with open(filepath1, 'r') as cov_report:
           header = cov_report.readline()
           header = header.rstrip()
           stats = cov_report.readline()
           stats = stats.rstrip()
           stats = stats.split()
           ref_name = stats[0]
           #print(ref_name)
           start = stats[1]
           end = stats[2]
           reads_mapped = stats[3]
           cov_bases = stats[4]
           cov = stats[5]
           depth = stats[6]
           baseq = stats[7]
           #print(reads_mapped)
           mapq = stats[8]
        
       #Get number of raw reads
       proc_1 = subprocess.run('zcat ' + "${params.output}"+"/dengue4/"+"${samplename}"+"/"+"${samplename}" + '_1_humanclean.fastq.gz | wc -l', shell=True, capture_output=True, text=True, check=True)
       wc_out_1 = proc_1.stdout.rstrip()
       reads_1 = int(wc_out_1) / 4
       proc_2 = subprocess.run('zcat ' + "${params.output}"+"/dengue4/"+"${samplename}"+"/"+"${samplename}" + '_2_humanclean.fastq.gz | wc -l', shell=True, capture_output=True, text=True, check=True)
       wc_out_2 = proc_2.stdout.rstrip()
       reads_2 = int(wc_out_2) / 4
       raw_reads = reads_1 + reads_2
       raw_reads = int(raw_reads)

       #Get number of clean reads
       proc_c1x = subprocess.run('zcat ' + "${params.output}"+"/dengue4/"+"${samplename}"+"/"+"${samplename}" + '_1.fq.gz | wc -l', shell=True, capture_output=True, text=True, check=True)
       wc_out_c1x = proc_c1x.stdout.rstrip()
       reads_c1x = int(wc_out_c1x) / 4
       proc_c2x = subprocess.run('zcat ' + "${params.output}"+"/dengue4/"+"${samplename}"+"/"+"${samplename}" + '_2.fq.gz | wc -l', shell=True, capture_output=True, text=True, check=True)
       wc_out_c2x = proc_c2x.stdout.rstrip()
       reads_c2x = int(wc_out_c2x) / 4
       clean_reads = reads_c1x + reads_c2x
       clean_reads = int(clean_reads)
       #print(clean_reads)
    
       #Get percentage of mapped reads/clean reads
       percent_map = "%0.4f"%((int(reads_mapped)/int(clean_reads))*100)
       #print(percent_map)
    
       #Gather QC metrics for consensus assembly
       filepath2 = "${params.output}"+"/dengue4/"+"${samplename}"+"/assembly/"+"${samplename}"+".consensus.fa"
       with open(filepath2, 'r') as assem:
           header = assem.readline()
           header = header.rstrip()
           bases = assem.readline()
           bases = bases.rstrip()
           num_bases = len(bases)
           ns = bases.count('N')
           called = num_bases - ns
           pg = "%0.4f"%((called/int(end))*100)
           #print(called)
           #print(end)
           #print(pg)
       #Rename header in fasta to just sample name
       subprocess.run("sed -i \'s/^>.*/>"+"${samplename}"+"/\' "+filepath2, shell=True, check=True)
    
       #QC flag
       pg_flag = ''
       dp_flag = ''
       qc_flag = ''
       if float(pg) < 79.5:
           pg_flag = 'FAIL: Percent genome < 80%'
           qc_flag = qc_flag + pg_flag
       else:
           if float(depth) < 100:
              dp_flag = 'FAIL: Mean read depth < 100x'
              qc_flag = qc_flag + dp_flag
           if qc_flag == '':
              qc_flag = qc_flag + 'PASS'
       #print(qc_flag)
       
       idname = "${samplename}".lstrip("SER4_")
       if qc_flag == 'PASS':      
           subprocess.run("cp "+filepath2+" "+"${params.output}"+"/dengue4/assemblies/"+idname+".consensus.fa", shell=True, check=True)   
           subprocess.run('cp ' + "${params.output}"+"/dengue4/"+"${samplename}" + '/variants/'+"${samplename}" + '.variants.tsv ' + "${params.output}"+'/dengue4/variants/'+idname+'.variants.tsv', shell=True, check=True)
    
           #Run VADR
           out_log = open("${params.output}"+"/dengue4/"+"${samplename}"+"/"+"${samplename}"+'.out', 'w')
           err_log = open("${params.output}"+"/dengue4/"+"${samplename}"+"/"+"${samplename}"+'.err', 'w')
           subprocess.run("singularity exec -B "+"${params.output}"+"/dengue4/"+"${samplename}"+"/assembly:/data docker://staphb/vadr:1.5.1 /opt/vadr/vadr/miniscripts/fasta-trim-terminal-ambigs.pl --minlen 50 --maxlen 30000 /data/" + "${samplename}"+".consensus.fa > " + "${params.output}"+"/dengue4/"+"${samplename}"+"/assembly/"+"${samplename}"+".trimmed.fasta", shell=True, stdout=out_log, stderr=err_log, check=True)

           subprocess.run("singularity exec -B "+"${params.output}"+"/dengue4/"+"${samplename}"+"/assembly:/data docker://staphb/vadr:1.5.1 v-annotate.pl --split --cpu 8 --group Dengue --nomisc --noprotid --mkey flavi --mdir /opt/vadr/vadr-models-flavi/ --noseqnamemax /data/"+"${samplename}"+".trimmed.fasta /data/"+"vadr_results", shell=True, stdout=out_log, stderr=err_log, check=True)


           #Parse through VADR outputs to get PASS or REVIEW flag
           vadr_flag = ''
           with open("${params.output}"+"/dengue4/"+"${samplename}"+"/assembly/vadr_results/vadr_results.vadr.pass.list", 'r') as p_list:
               result = p_list.readline()
               result = result.rstrip()
               if result == "${samplename}":
                   vadr_flag = 'PASS'

           with open("${params.output}"+"/dengue4/"+"${samplename}"+"/assembly/vadr_results/vadr_results.vadr.fail.list", 'r') as f_list:
               result = f_list.readline()
               result = result.rstrip()
               if result == "${samplename}":
                   vadr_flag = 'REVIEW'

           #Copy VADR error report to main analysis folder for easier review
           if vadr_flag == 'REVIEW':
               subprocess.run("cp " + "${params.output}"+"/dengue4/"+"${samplename}"+"/assembly/vadr_results/vadr_results.vadr.alt.list "+"${params.output}"+"/dengue4/vadr_error_reports/"+"${samplename}"+".vadr.alt.list", shell=True, check=True)

           out_log.close()
           err_log.close()
           subprocess.run("mv " + "${params.output}"+"/dengue4/"+"${samplename}"+"/assembly/vadr_results "+"${params.output}"+"/dengue4/"+"${samplename}"+"/assembly/"+"${samplename}"+"_vadr_results", shell=True, check=True)

       else:
           vadr_flag = 'NA'
           lineage = 'NA'
           sotc_out = 'NA'
        
       with open("${params.output}"+"/dengue4/"+"${samplename}"+"/report.txt", 'w') as report:
           header = ['sampleID', 'reference', 'start', 'end', 'num_raw_reads', 'num_clean_reads', 'num_mapped_reads', 'percent_mapped_clean_reads', 'cov_bases_mapped', 'percent_genome_cov_map', 'mean_depth', 'mean_base_qual', 'mean_map_qual', 'assembly_length', 'numN', 'percent_ref_genome_cov', 'VADR_flag', 'QC_flag']
           report.write('\t'.join(map(str,header)) + '\n')
           results = [idname, ref_name, start, end, raw_reads, clean_reads, reads_mapped, percent_map, cov_bases, cov, depth, baseq, mapq, num_bases, ns, pg, vadr_flag, qc_flag]
           report.write('\t'.join(map(str,results)) + '\n')

    else:
       print("check pystats.nf step.")
    
    /$
}