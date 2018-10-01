# v3seq: quick reconstruction of V3 haplotypes

`v3seq` is a quick pipeline that extracts reads overlapping V3 by aligning
them with blastx-fast to a set of V3 consensus sequences (taken from
[this paper](https://doi.org/10.1016/S0165-2478(02)00101-3)). Reads are then
oriented in the correct direction, an MSA is produced and a rough estimate
of haplotypes is returned by counting identical sequences over the V3 region.
