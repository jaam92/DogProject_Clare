library(tidyverse)

#read files in

clare = read_delim("clareSNPs.bim", delim = "\t", col_names = c("chrom", "id", "pos"))
jaz = read_delim("jazSNPs.bim", delim = "\t", col_names = c("chrom", "id", "pos"))

indivs_jaz = read_delim("MergedFile_CornellCanineFitak_UnrelatedsOnly.fam", delim = " ", col_names = c("id", "breed"), col_types = cols_only(id = 'c', breed = 'c')) %>%
  group_by(breed) %>%
  top_n(20)

write.table(indivs_jaz, "extractUnrelatedsJazData_n20.txt", sep = "\t", row.names = F, col.names = F, quote = F)

#intersect them and add an id col
intersected = jaz %>%
  inner_join(clare, by=c("chrom", "pos"), suffix = c("_j", "_c"))


keep = intersected %>%
  select(id_c)

write.table(keep, "extractSNPs.txt", sep = "\t", row.names = F, col.names = F, quote = F)

renameSNPs = intersected %>%
  select(id_c, id_j)

write.table(renameSNPs, "updateRSID.txt", sep = "\t", row.names = F, col.names = F, quote = F)


