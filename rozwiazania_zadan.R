## Zad1

library(stringr)

upperlower <- function (x) {
  x_split <- unlist(strsplit(x, ''))
  even_list <- list()
  odd_list <- list()
  for (i in 1:length(x_split)) {
    if (i %% 2 == 0) {
      even_list <- append(even_list, toupper(x_split[i]))
    } else {
      even_list <- append(even_list, x_split[i])
    } 
  }
  for (i in 1:length(x_split)) {
    if (i %% 2 != 0) {
      odd_list <- append(odd_list, toupper(x_split[i]))
    } else {
      odd_list <- append(odd_list, x_split[i])
    } 
  }
return(list(
  paste(even_list, collapse = ''),
  paste(odd_list,collapse = '')
))
}

test <- upperlower(letters[5:20])
test

## Zad2

test1 <- 'ABBA'
test2 <- 'aBcbA'
test3 <- 'RabarbArka'

letter_count <- function (x) {
  x_split <- unlist(strsplit(x, ''))
  low <- tolower(x_split)
  t <- table(low)
  for (i in 1:dim(t)) {
    if (t[[i]] > 1) {
      print(
        t[i]
        )
      }
  }
}

letter_count(test3)

### Zad 3 (mo?nna to szybciej zrobi? tabixem ale rozumiem ?e sprawdzian zak?ada korzystanie z R)

library(data.table)
vcf <- fread('../CPCT02220079.annotated.processed.vcf.gz',skip=402,sep='\t')

colnames(vcf)[1] <- 'Chromosome'

library(tidyverse)

chr12 <-  vcf %>% filter(Chromosome == '12' & POS >= 112204691 & POS <= 112247789)
write.table(chr12,'chr12_010421.csv',quote = F,row.names = F,sep='\t')

### Zad 4

### Length > 0 = insercje
### Length < 0 = delecje
vcf$Length <- nchar(vcf$REF) - nchar(vcf$ALT)

plot_vcf <- vcf %>% group_by(Chromosome,Length) %>% summarise(Count = log(n())) %>% filter(Length != 0 & Count > 0) %>% 
  mutate(Type = ifelse(Length < 0,'Insertion','Deletion'))

plot_vcf  %>% ggplot() + geom_col(aes(x=Length, y= Count, fill=Type)) + 
  facet_wrap(~Chromosome) + ggsave('indels.png',height = 15, width=15) + ylab('log10 count')

write.table(plot_vcf,'indels_010421.csv',quote = F,row.names = F,sep='\t')

## Zad5
vcf_het <- vcf %>% filter(FILTER == 'PASS')
vcf_het$GoNLv5 <- str_extract(vcf_het$INFO, ';GoNLv5_AF=.*;')
vcf_het$GoNLv5 <- gsub('GoNLv5_AF=','',vcf_het$GoNLv5)
vcf_het$GoNLv5 <- as.numeric(gsub(';','',vcf_het$GoNLv5))

vcf_het$af.05 <- str_extract(vcf_het$INFO, ';AF=0.5;')
vcf_het$af.05 <- gsub(';AF=','',vcf_het$af.05)
vcf_het$af.05 <- as.numeric(gsub(';','',vcf_het$af.05))

vcf_het_filtered <- vcf_het %>% filter(is.na(af.05) == F) %>% filter(GoNLv5 < 0.01 & is.na(GoNLv5) == F)

### Liczba wariant?w z AF < 0.01
vcf_het_filtered %>% summarise('AF < 0.01'=n())


## Zad6

noauto <- c('X','Y','MT')

avg_cov <- vcf %>% filter(!(Chromosome %in% noauto)) %>% 
  separate(CPCT02220079R, c('GT', 'AD', 'DP'), sep=':',extra='drop') %>% select(Chromosome, DP)

avg_cov$DP <- as.numeric(avg_cov$DP)

avg_cov[is.na(avg_cov$DP)] <- 0
mean(avg_cov$DP)

chr_plot <- avg_cov %>% group_by(Chromosome) %>% summarise(avg_cov = mean(DP))
chr_plot$Chromosome <- as.numeric(chr_plot$Chromosome)

chr_plot %>% ggplot(aes(Chromosome,avg_cov)) + geom_bar(stat = 'identity') + ylab('Mean depth of coverage') + 
  xlab('Chromosomes') + scale_x_continuous(breaks=seq(1,22,1)) + ggsave('cov_per_chr.png',height = 15, width=15)  

write.table(chr_plot, 'cov_per_chr.csv',sep='\t',quote = F,row.names = F)

## Zad 8

### Funkcja do obliczania VAF dla indeli
library(vcfR)
strelka_vaf <- function(x) {
  streika_vcf <- read.vcfR(x)
  streika_tidy <- vcfR2tidy(streika_vcf)
  streika_indel <- streika_tidy$gt
  
  streika_vaf <- streika_indel %>% select(Indiv, gt_TAR,gt_TIR)
  streika_vaf <- streika_vaf %>% separate(gt_TAR, into = c('tier1RefCounts','tar'),remove = T,convert = T) %>% select(-tar)
  streika_vaf <-  streika_vaf %>% separate(gt_TIR, into = c('tier1AltCounts','tir'),remove = T,convert = T) %>% select(-tir)
  
  streika_vaf$VAF <- streika_vaf$tier1AltCounts/(streika_vaf$tier1RefCounts+streika_vaf$tier1AltCounts)
  streika_vaf$VAF <- ifelse(is.nan(streika_vaf$VAF),0,streika_vaf$VAF)
  return(streika_vaf)
}

### VAF dla indeli
indel <- strelka_vaf('../T1_vs_N1_head.strelka.somatic.indels.norm.vcf.gz')
write.table(indel,'streika_indel.csv',sep='\t',quote=F,row.names = F)

### VAF dla snv

#### Nie wiedzia?em jak to zrobi? (znaczy wiedziałem ale by??o to zbyt toporne by udost?spnia? to publicznie) 
#### wi?c postanowi?em wykorzysta? gotowe rozwi?zanie https://github.com/biobenkj/StrelkaParser

parse_vcf_alt1 <- function(x) {
  vcf <- readr::read_tsv(x, comment = "##")
  vcf <- dplyr::rename(vcf, CHROM = `#CHROM`)
  vcf <- dplyr::mutate(vcf, ALT = strsplit(ALT, ",") %>% purrr::map_chr(1),
                       info = parse_info(INFO),
                       normal = parse_format(FORMAT, NORMAL),
                       tumour = parse_format(FORMAT, TUMOR))
  vcf <- tidyr::unnest(vcf, info)
  vcf <- tidyr::unnest(vcf, N = normal, T = tumour, .sep = "_")
  vcf <- purrr::lmap_at(vcf, dplyr::matches("[TN]_[ACTG]U", vars = names(vcf)), split_tiers)
  vcf <- readr::type_convert(vcf)
  vcf <- tidyr::gather(vcf, allele, count, dplyr::matches("[TN]_[ACTG]U_[12]"))
  vcf <- tidyr::separate(vcf, allele, c("sample", "base", "tier"), sep = "[_U]+", remove = FALSE)
  vcf <- dplyr::group_by(vcf, CHROM, POS, REF, ALT)
  vcf <- dplyr::mutate(vcf, T_REF_COUNT = count[REF == base & sample == "T" & tier == TQSS_NT],
                       T_ALT_COUNT = count[ALT == base & sample == "T" & tier == TQSS_NT],
                       N_REF_COUNT = count[REF == base & sample == "N" & tier == TQSS_NT],
                       N_ALT_COUNT = count[ALT == base & sample == "N" & tier == TQSS_NT],
                       T_VAF = T_ALT_COUNT / (T_ALT_COUNT + T_REF_COUNT),
                       N_VAF = N_ALT_COUNT / (N_ALT_COUNT + N_REF_COUNT))
  vcf <- dplyr::select(vcf, -sample, -base, -tier)
  vcf <- tidyr::spread(vcf, allele, count)  # Takes up most of the time
  vcf
}
row_tibble <- function(x, col_names) {
  tibble::as_tibble(rbind(setNames(x, col_names)))
}

parse_format <- function(format, sample) {
  purrr::map2(
    strsplit(sample, ":", fixed = TRUE),
    strsplit(format, ":", fixed = TRUE),
    row_tibble)
}

parse_info <- function(info) {
  strsplit(info, ";", fixed = TRUE) %>%
    purrr::map(~row_tibble(sub("^.*=(.*)", "\\1", .x), sub("^(.*)=.*", "\\1", .x)))
}

split_tiers <- function(x) {
  name <- names(x)
  x$`1` <- x[[1]] %>%
    strsplit(",", fixed = TRUE) %>%
    purrr::map(~as.integer(.x[1])) %>%
    purrr::flatten_int()
  x$`2` <- x[[1]] %>%
    strsplit(",", fixed = TRUE) %>%
    purrr::map(~as.integer(.x[2])) %>%
    purrr::flatten_int()
  names(x)[2:3] <- paste0(name, "_", names(x)[2:3])
  x[2:3]
}

calc_counts <- function(x) {
  T_REF_COUNT_KEY <- paste0("T_", x$REF, "U_", x$TQSS_NT)
  T_ALT_COUNT_KEY <- paste0("T_", x$ALT, "U_", x$TQSS_NT)
  N_REF_COUNT_KEY <- paste0("N_", x$REF, "U_", x$TQSS_NT)
  N_ALT_COUNT_KEY <- paste0("N_", x$ALT, "U_", x$TQSS_NT)
  T_REF_COUNT <- x[[T_REF_COUNT_KEY]]
  T_ALT_COUNT <- x[[T_ALT_COUNT_KEY]]
  N_REF_COUNT <- x[[N_REF_COUNT_KEY]]
  N_ALT_COUNT <- x[[N_ALT_COUNT_KEY]]
  T_VAF <- T_ALT_COUNT / (T_ALT_COUNT + T_REF_COUNT)
  N_VAF <- N_ALT_COUNT / (N_ALT_COUNT + N_REF_COUNT)
  data.frame(T_VAF = T_VAF, N_VAF = N_VAF)
}

snv <- parse_vcf_alt1('../T1_vs_N1_head.strelka.somatic.snvs.norm.vcf.gz')
snv_save <- snv %>% select(CHROM,POS,ID,REF,ALT,T_VAF,N_VAF)
write.table(snv_save,'streika_snv.csv',sep='\t',quote = F,row.names = F)


