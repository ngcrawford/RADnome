### Welcome to RADnome.org

Initial computational methods for analyzing restriction site associated DNA markers ([RAD-seq][1]) were designed for either single end sequencing or for contigs assembled from single digest RADs. However, new methods using double digests (e.g., [ddRAD-seq][2]) produce RADs where individual PE reads are non-overlapping.

To address this problem this software independently infers contigs for each R1 and R2 read using [rainbow][3]. Contigs representing read pairs are then associated based on read mapping positions. After contigs are paired all contigs are merged into a pseudo-genome (= RADnome). Individual samples are then mapped to this RADnome and software packages such as [GATK][4] can be used to identify SNVs.

### Schedule

RADnome is still in development. However, we expect to have a working pipeline in place by July 1st 2013 at the latest. Follow us on twitter or github and we'll let you know when everything is ready to go. 

### Contact

*  Email: RADnome@gmail.com
*  Twitter: [RADnome][5]
*  Github: [github.com/ngcrawford/RADnome][6]


### Authors and Contributors
This project was conceived W. Brian Simison (@wbsimey) and Nicholas Crawford (@ngcrawford) with contributions from Josh Pollock (@iliketurkey).


[1]: http://en.wikipedia.org/wiki/Restriction_site_associated_DNA_markers
[2]: http://www.plosone.org/article/info:doi/10.1371/journal.pone.0037135
[3]: dx.doi.org/10.1093/bioinformatics/bts482
[4]: http://www.broadinstitute.org/gatk/
[5]: https://twitter.com/RADnome
[6]: https://github.com/ngcrawford/RADnome
