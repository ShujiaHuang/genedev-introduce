# 基因组学竞赛题目

[![logo](https://github.com/gcta/genedev-introduce/blob/master/assets/GCTA.logo.png)](http://www.gene-dev.cn/)

该竞赛的相关内容已同步发布于GCTA官网[链接](http://www.gene-dev.cn/navigation/show_navigation?nid=17060)

## 1 数据存储问题：FASTQ数据的高倍数压缩

**【背景】**随着测序技术的进步，全世界测序产出的DNA数据正在激增中，如何有效降低[DNA](https://zh.wikipedia.org/wiki/%E8%84%B1%E6%B0%A7%E6%A0%B8%E7%B3%96%E6%A0%B8%E9%85%B8)数据的存储空间已经成为了一个急需解决的难题。目前一般的压缩格式（gz，bz2等）只能将DNA的测序数据（FASTQ格式，或称fq格式）压缩至原来的30%左右。但fq数据有着自己固定的格式形式，DNA测序数据也只有4种碱基，由A，C，G，T这四个字母表示，是能够针对其数据上的特点实现更高倍数的压缩的。

**【题目】**给定一个[FASTQ格式](https://zh.wikipedia.org/wiki/FASTQ%E6%A0%BC%E5%BC%8F)的DNA测序数据文件，非压缩状态下该文件的大小约为150GB，现要求将其至少无损压缩至原来的1/15。本题所需数据[下载地址]()。

【要求】

- 原创或改进现有方法；
- 无损压缩，信息不能有任何丢失；
- 通用性，要求对任意的FASTQ数据有效；
- 时效性，压缩和解压时间必须在可接受范围内，不能超过传统压缩方法（[gzip](https://zh.wikipedia.org/wiki/Gzip)）的3倍。

**【评核标准】:**

- 压缩比（60%）；
- 压缩和解压时间（40%）。

**【加分项】**

- 其它有助于读写该压缩格式的附加功能给予适当加分（10）

**【FASTQ文件的格式说明】**FASTQ是当前存储物种的原始测序数据（DNA）的标准文件，它每四行为一个独立的序列存储单元（成为read单元），每一个read单元的格式如下：

- 序列标识以及相关的描述信息，以‘@’开头，是每一个序列的唯一标识符；
- 第二行是序列，由A，C，G，T和N构成，以当前高通量测序技术其长度一般在100个字符以上，不同测序仪会有差别，其中A，C，G，T是碱基信息，N为测序失败时用来替补的补位码；
- 第三行以‘+’开头，后面是序列标示符、描述信息，或者什么也不加（本题所用的数据该行只有「+」，评估也将以只有「+」的这类情况评定）
- 第四行，是[序列的质量信息](https://en.wikipedia.org/wiki/Phred_quality_score)，和第二行序列中的碱基一一对应，每一个碱基对应一个质量值，质量值用ASCII码表示，用以衡量该测序碱基的可靠程度，质量值越高越可靠。


例如：

```
@ERR194147.1 HSQ1004:134:C0D8DACXX:1:1104:3874:86238/1
GGTTCCTACTTCAGGGTCATAAAGCCTAAATAGCCCACACGTTCCCCTTAAATAAGACATCACGATGGATCACAGGTCTATCACCCTATTAACCACTCACG
+
CC@FFFFFHHHHHJJJFHIIJJJJJJIHJIIJJJJJJJJIIGIJJIJJJIJJJIJIJJJJJJJJJJIJHHHHFFFDEEEEEEEEDDDCDDEEDDDDDDDDD
```


## 2 计算问题：高效的人类全基因数据分析

**【背景】**人的[基因组](https://zh.wikipedia.org/wiki/%E4%BA%BA%E9%A1%9E%E5%9F%BA%E5%9B%A0%E7%B5%84)为3G，用于[全基因组测序](https://en.wikipedia.org/wiki/Whole_genome_sequencing)数据分析时，需要[测序深度](http://baike.baidu.com/link?url=9QQtc999YINr7u5ExQ-YPWn3SoktRuGPNYQnZ9m3luaUgASenKuVrLjCZuBg_x7404i4pPMxghR8fVjINbhkUq)为50x或者更高，使用常见的生物信息分析工具和方法（[bwa](http://bio-bwa.sourceforge.net/)+[picard](http://www.psc.edu/index.php/user-resources/software/picard)+[GATK](https://www.broadinstitute.org/gatk/)），时间基本需要10个-20个小时，这对于日益增长的人类基因组数据来说是远远不能满足数据解读的速度需求的。当前基因组数据分析的最大瓶颈是，数据的解读速度远不及数据的产出速度，全基因组数据分析是人类基因组数据解读中最基本的一个步骤。

**【题目】**30分钟完成50x-60x人类全基因组数据标准分析（从fq数据到变异数据的产出），本题提供的数据约覆盖人类基因组55x。本题所需数据[下载地址]()。

**【要求】**

- 原创或改进现有方法；
- 方法不限，但所用资源和成本应具有实际的可行性和可推广性；
- 方法和方案必须完整可复用；
- 必须是Pair-End测序的数据，如本题所提供；
- 整个分析流程不局限于目前全基因组数据分析的传统过程[1], 但最后必须要以标准的VCF/gVCF格式输出变异数据；
- 本题所用的人类参考序列版本统一为：GRCh37（或称hg19），为了确保参赛者们所使用参考序列的一致性，组委会统一提供下载[链接](ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.gz)；
- 流程监控，可以借助外部工具，或者流程内置相应的监控工具等方式，最基本的要求是，必须能够准确监控流程中每个步骤的运行状态（成功/失败），并返回相应的处理值，以便后续处理。

> [1] 传统的全基因组分析过程一般包括：原始fq数据质控和低质量数据的过滤，比对，碱基质量值重校正，变异区域重比对，SNP与Indel的变异检测，变异数据重校正等这些处理过程。

**【评核标准】**

- 时间长短，时间为用户完成分析的时间，非电脑计算所需的时间，比如竞赛者同时使用多台电脑进行并行分析，时间上并不累加这些并行运算的时间（30）；
- 资源和成本可接受性（20）；
- 方法的易用性和复用性（20）；
- 变异数据检测的准确性，与本例提供的变异集合进行比较（30）：
1）SNP的一致率（位置，碱基和Genotype同时一致）至少要到达99.5%；
2）Indel的类型，genotype和断点同时一致的比率至少98.5%；
3）SNP和Indel的假阴和假阳率都必须小于或等于1%。
- 流程监控的有效性，该项只要达到了最基本的要求便合格，没有该项功能，总分将直接扣10分；


> 对于该题的变异准确性评估，若参赛者能以充分的理由说明自己所检测出来的那些不一致或者假阳的SNP和Indel确实要比组委会所提供的更准确则可以不依上述的准确性评估标准！

**【本题所用数据说明】**
本题所提供的原始测序数据（FASTQ格式）为Pair-End测序数据，来自[千人基因组计划](http://www.1000genomes.org/)中的一个样本（NA12878），本次竞赛分析所需的基本数据只有3个，具体如下：

```
1. NA12878_1.fastq.gz  NA12878的read1，长度101bp，碱基质量值为ASCII-33体系
2. NA12878_1.fastq.gz  NA12878的read2，长度101bp，碱基质量值为ASCII-33体系
3. human_g1k_v37.fasta.gz 人类基因组参考序列
```




