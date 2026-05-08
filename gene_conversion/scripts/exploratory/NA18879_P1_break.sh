#  Sample:          NA18879_
#  Assembly:        data/masked/masked_sequences/NA18879_masked.fa
#  Query:           work/P1Long/arms/HG00621_A.fa
#  Preset:          asm5
#  qcov min:        20.0%
#  Indel min:       5 bp
#  Outdir:          output_missing_P1
#  
#  Raw PAF saved: output_missing_P1/NA18879__detailed.paf
#  
#  Detailed alignment report — NA18879_
#  Query: work/P1Long/arms/HG00621_A.fa  (1,204,151 bp)
#  Assembly: data/masked/masked_sequences/NA18879_masked.fa
#  ────────────────────────────────────────────────────────────────────────────────
#  #    contig                         strand     t_start      t_end    qcov  pident  mapq  indels
#  ────────────────────────────────────────────────────────────────────────────────
#  1    NA18879_chrY                   minus   25,733,045 26,937,202  100.0%   99.9%    60      20
#                                          deletion      9 bp
#                                         insertion      5 bp
#                                         insertion      8 bp
#                                          deletion      7 bp
#                                          deletion     10 bp
#                                          deletion      8 bp
#                                          deletion      9 bp
#                                          deletion     10 bp
#                                          deletion     10 bp
#                                         insertion      6 bp
#                                         insertion      6 bp
#                                          deletion     14 bp
#                                         insertion     10 bp
#                                         insertion      5 bp
#                                         insertion      6 bp
#                                         insertion     14 bp
#                                         insertion     12 bp
#                                          deletion      6 bp
#                                         insertion     11 bp
#                                          deletion      8 bp
#  2    NA18879_chrY                   plus    24,986,665 25,633,984   53.8%   99.9%     0      15
#                                          deletion      7 bp
#                                         insertion     15 bp
#                                         insertion      6 bp
#                                         insertion      6 bp
#                                          deletion      8 bp
#                                          deletion      8 bp
#                                          deletion      7 bp
#                                         insertion     12 bp
#                                          deletion     12 bp
#                                         insertion     18 bp
#                                          deletion      8 bp
#                                         insertion      5 bp
#                                         insertion    101 bp
#                                          deletion      9 bp
#                                         insertion      8 bp
#  3    NA18879_chrY                   plus    24,432,374 24,982,720   45.7%   99.9%     0      13
#                                         insertion      8 bp
#                                         insertion      7 bp
#                                         insertion     11 bp
#                                         insertion     12 bp
#                                          deletion      8 bp
#                                          deletion      8 bp
#                                         insertion      6 bp
#                                          deletion      8 bp
#                                         insertion      8 bp
#                                          deletion     10 bp
#                                         insertion      6 bp
#                                         insertion      6 bp
#                                          deletion      8 bp
#  4    NA18879_chrY                   plus    23,636,353 24,007,659   30.8%   99.9%     0      13
#                                         insertion      6 bp
#                                         insertion     11 bp
#                                          deletion     20 bp
#                                         insertion      6 bp
#                                         insertion     28 bp
#                                          deletion     12 bp
#                                         insertion      5 bp
#                                          deletion     12 bp
#                                         insertion      6 bp
#                                          deletion     13 bp
#                                          deletion     11 bp
#                                          deletion      9 bp
#                                         insertion     10 bp
#  ────────────────────────────────────────────────────────────────────────────────
#  Total hits reported: 4
#  Flags (standard thresholds): hit_count=1
#  
#  Report saved: output_missing_P1/NA18879__detailed_report.txt



samtools faidx /Users/kxp5629/proj/HPRC_Y/data/masked/masked_sequences/NA18879_masked.fa NA18879_chrY:25733045-26937202 > NA18879_minus.fa
samtools faidx /Users/kxp5629/proj/HPRC_Y/data/masked/masked_sequences/NA18879_masked.fa NA18879_chrY:24432374-25633984 > NA18879_plus.fa