Diversity\_Analysis
================
Nick Allsing
3/3/2022

\#\#\#Import OTU table into R for analysis I am using R version 4.0.4
(2021-02-15) Now the newly created and modified metagenome (or bacteria,
archaea, eukarya, virus) file can be imported into R for analysis

``` r
library(propr)
library(magrittr)
library(zCompositions)
```

    ## Loading required package: MASS

    ## Loading required package: NADA

    ## Loading required package: survival

    ## 
    ## Attaching package: 'NADA'

    ## The following object is masked from 'package:stats':
    ## 
    ##     cor

    ## Loading required package: truncnorm

``` r
library(dplyr)
```

    ## 
    ## Attaching package: 'dplyr'

    ## The following object is masked from 'package:MASS':
    ## 
    ##     select

    ## The following objects are masked from 'package:stats':
    ## 
    ##     filter, lag

    ## The following objects are masked from 'package:base':
    ## 
    ##     intersect, setdiff, setequal, union

``` r
library(ggplot2)
library(tidyverse) 
```

    ## ── Attaching packages ─────────────────────────────────────── tidyverse 1.3.1 ──

    ## ✓ tibble  3.1.0     ✓ purrr   0.3.4
    ## ✓ tidyr   1.1.3     ✓ stringr 1.4.0
    ## ✓ readr   1.4.0     ✓ forcats 0.5.1

    ## ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
    ## x tidyr::extract()   masks magrittr::extract()
    ## x dplyr::filter()    masks stats::filter()
    ## x dplyr::lag()       masks stats::lag()
    ## x dplyr::select()    masks MASS::select()
    ## x purrr::set_names() masks magrittr::set_names()
    ## x purrr::simplify()  masks propr::simplify()

``` r
library(readr)
library(mixOmics)
```

    ## Warning: package 'mixOmics' was built under R version 4.0.5

    ## Loading required package: lattice

    ## 
    ## Attaching package: 'lattice'

    ## The following object is masked from 'package:propr':
    ## 
    ##     parallel

    ## 
    ## Loaded mixOmics 6.14.1
    ## Thank you for using mixOmics!
    ## Tutorials: http://mixomics.org
    ## Bookdown vignette: https://mixomicsteam.github.io/Bookdown
    ## Questions, issues: Follow the prompts at http://mixomics.org/contact-us
    ## Cite us:  citation('mixOmics')

    ## 
    ## Attaching package: 'mixOmics'

    ## The following object is masked from 'package:purrr':
    ## 
    ##     map

    ## The following object is masked from 'package:propr':
    ## 
    ##     pca

``` r
library(selbal)
library(vegan)
```

    ## Loading required package: permute

    ## This is vegan 2.5-7

``` r
library(scales)
```

    ## 
    ## Attaching package: 'scales'

    ## The following object is masked from 'package:purrr':
    ## 
    ##     discard

    ## The following object is masked from 'package:readr':
    ## 
    ##     col_factor

``` r
library(ggthemes)
require(gtools)
```

    ## Loading required package: gtools

    ## 
    ## Attaching package: 'gtools'

    ## The following object is masked from 'package:permute':
    ## 
    ##     permute

``` r
library(viridis)
```

    ## Loading required package: viridisLite

    ## 
    ## Attaching package: 'viridis'

    ## The following object is masked from 'package:scales':
    ## 
    ##     viridis_pal

``` r
library(ggpubr)
library(corrplot)
```

    ## corrplot 0.84 loaded

``` r
library(psych)
```

    ## 
    ## Attaching package: 'psych'

    ## The following object is masked from 'package:gtools':
    ## 
    ##     logit

    ## The following objects are masked from 'package:scales':
    ## 
    ##     alpha, rescale

    ## The following object is masked from 'package:mixOmics':
    ## 
    ##     pca

    ## The following objects are masked from 'package:ggplot2':
    ## 
    ##     %+%, alpha

    ## The following object is masked from 'package:propr':
    ## 
    ##     pca

``` r
library(rstatix)
```

    ## 
    ## Attaching package: 'rstatix'

    ## The following object is masked from 'package:MASS':
    ## 
    ##     select

    ## The following object is masked from 'package:stats':
    ## 
    ##     filter

``` r
library(png)
library(grid)
library(gridExtra)
```

    ## 
    ## Attaching package: 'gridExtra'

    ## The following object is masked from 'package:dplyr':
    ## 
    ##     combine

``` r
library(cowplot)
```

    ## 
    ## Attaching package: 'cowplot'

    ## The following object is masked from 'package:ggpubr':
    ## 
    ##     get_legend

    ## The following object is masked from 'package:ggthemes':
    ## 
    ##     theme_map

``` r
library(data.table)
```

    ## 
    ## Attaching package: 'data.table'

    ## The following object is masked from 'package:purrr':
    ## 
    ##     transpose

    ## The following objects are masked from 'package:dplyr':
    ## 
    ##     between, first, last

``` r
# read in OTU table and mapping file
metagen_OTU <- read.csv('./metagenomemodified.csv',sep=',',row.names = 1,header = 1,check.names=FALSE)
map <- read.table('./TJRIVER_mappingfile.txt',sep='\t',header =1,row.names = 1)
```

``` r
#To reduce the overall number of variables if desired.
metagen_1 <- metagen_OTU[rowSums(metagen_OTU <1) <= 0.01 * ncol(metagen_OTU), ]
metagen_5 <- metagen_OTU[rowSums(metagen_OTU <1) <= 0.05 * ncol(metagen_OTU), ]
metagen_8 <- metagen_OTU[rowSums(metagen_OTU <1) <= 0.08 * ncol(metagen_OTU), ]
metagen_10 <- metagen_OTU[rowSums(metagen_OTU <1) <= 0.10 * ncol(metagen_OTU), ]
metagen_50 <- metagen_OTU[rowSums(metagen_OTU <1) <= 0.50 * ncol(metagen_OTU), ]
metagen_100 <- metagen_OTU[rowSums(metagen_OTU <1) <= 1.00 * ncol(metagen_OTU), ]
```

``` r
#Check for number of 0s in sample
table(metagen_OTU == 0)
```

    ## 
    ##  FALSE   TRUE 
    ## 403193 218021

``` r
#replace 0s with p-counts
metagen_z <- cmultRepl((metagen_OTU), output = 'p-counts')
```

    ## No. corrected values:  218017

``` r
#Check to see if 0s have been removed
table(metagen_z == 0)
```

    ## 
    ##  FALSE 
    ## 621214

``` r
#Transpose for clr
metagen_z_t <- t(metagen_z)
```

``` r
#define clr function
clr <- function(x) sweep(log(x), 1, rowMeans(log(x)), "-")
```

``` r
#Perform clr
metagen_clr_t <- clr(metagen_z_t)
```

``` r
#transpose back
metagen_clr <- t(metagen_clr_t)
```

``` r
#Create dataframes and order them and map by ID
metagen_clr_df <- as.data.frame(metagen_clr)
metagen_clr_ordered <- metagen_clr_df[order(row.names(metagen_clr_df)), order(colnames(metagen_clr_df))]
map_ordered <- map[order(row.names(map)),order(colnames(map))]
```

``` r
#Write out clr file
write.csv(metagen_clr_ordered, "./TJ_River_metagen_clr.csv")
```

``` r
#Filter out anything not on mapping file
metagen_clr_filtered <- metagen_clr_ordered[,row.names(map_ordered)]
```

``` r
#Transpose data
metagen_clr_filtered_t <- t(metagen_clr_filtered)
```

``` r
#Make matrix
metagen_clr_m <- as.matrix(metagen_clr_filtered_t)
```

``` r
#Make NMDS data
metagen_nmds = metaMDS(metagen_clr_m, distance = 'euclidean')
```

    ## 'comm' has negative data: 'autotransform', 'noshare' and 'wascores' set to FALSE

    ## Run 0 stress 0.1123695 
    ## Run 1 stress 0.116438 
    ## Run 2 stress 0.1248176 
    ## Run 3 stress 0.1165616 
    ## Run 4 stress 0.1316443 
    ## Run 5 stress 0.1227107 
    ## Run 6 stress 0.1246235 
    ## Run 7 stress 0.1226695 
    ## Run 8 stress 0.1236695 
    ## Run 9 stress 0.1180747 
    ## Run 10 stress 0.1169033 
    ## Run 11 stress 0.1190135 
    ## Run 12 stress 0.1142948 
    ## Run 13 stress 0.1224159 
    ## Run 14 stress 0.1218124 
    ## Run 15 stress 0.1150836 
    ## Run 16 stress 0.1228996 
    ## Run 17 stress 0.1161725 
    ## Run 18 stress 0.1286614 
    ## Run 19 stress 0.1250439 
    ## Run 20 stress 0.1179974 
    ## *** No convergence -- monoMDS stopping criteria:
    ##     19: stress ratio > sratmax
    ##      1: scale factor of the gradient < sfgrmin

``` r
#Extract NMDS data
metagen_scores = as.data.frame(scores(metagen_nmds))
```

``` r
#Add mapping file columns to NMDS data
metagen_scores$date = as.factor(map_ordered$date)
metagen_scores$number = as.factor(map_ordered$number)
metagen_scores$site = as.factor(map_ordered$site)
```

``` r
#plot NMDS adjust aes according to data
metagen_nmdsplot <- ggplot(metagen_scores, aes(x = NMDS1, y = NMDS2)) +
  geom_point(size = 3, aes(colour = site, shape = date))+ 
  geom_text(data = metagen_scores,aes(x=NMDS1,y=NMDS2,label=number),size=2,vjust=2)+
  ggtitle('Metagenomic NMDS')+
  theme(
    plot.title = element_text(color="black", size=10, face="bold"),
    panel.background = element_rect(fill = "white",
                                    colour = "black",
                                    size = 0.5, linetype = "solid"),
  )
```

``` r
#show plot
printplot <- ggpubr::ggarrange(metagen_nmdsplot,ncol = 1,nrow=1,labels = c("A"))
printplot
```

![](Diversity_Analysis_files/figure-gfm/unnamed-chunk-20-1.png)<!-- -->

``` r
ggsave('./NMDS_TJ_River_metagen.png', width=6.5, height=4,plot=metagen_nmdsplot)
```

``` r
#Vegdist
metag_dist<-vegdist(metagen_clr_filtered_t, method='euclidian')
```

``` r
#adonis
metag_site <- adonis2(metag_dist~site, data=map_ordered, permutations = 9999, method = "euclidian")
metag_date <- adonis2(metag_dist~date, data=map_ordered, permutations = 9999, method = "euclidian")
metag_number <- adonis2(metag_dist~number, data=map_ordered, permutations = 9999, method = "euclidian")
```

``` r
capture.output(metag_site,file="./TJ_River_site_permanova_10000.csv")
capture.output(metag_date,file="./TJ_River_date_permanova_10000.csv")
capture.output(metag_number,file="./TJ_River_number_permanova_10000.csv")
```

``` r
site_metagen <- factor(map_ordered[colnames(metagen_clr_filtered),]$site)
date_metagen <- factor(map_ordered[colnames(metagen_clr_filtered),]$date)
number_metagen <- factor(map_ordered[colnames(metagen_clr_filtered),]$number)
```

``` r
metagbeta_clr_site <- betadisper(metag_dist, site_metagen)
metagbetadispr_site <- permutest(metagbeta_clr_site)
metagdistances_site <- data.frame(metagbeta_clr_site[['distances']])
metagdistances_site <- cbind (metagdistances_site, site_metagen)
colnames(metagdistances_site) <- c('Distance_to_Centroid', 'site')
```

``` r
metagbeta_clr_date <- betadisper(metag_dist, date_metagen)
metagbetadispr_date <- permutest(metagbeta_clr_date)
metagdistances_date <- data.frame(metagbeta_clr_date[['distances']])
metagdistances_date <- cbind (metagdistances_date, date_metagen)
colnames(metagdistances_date) <- c('Distance_to_Centroid', 'date')
```

``` r
dist_microbe_box_site <- ggboxplot(metagdistances_site, x = "site", y = "Distance_to_Centroid" , add = "point")
dist_microbe_box_site
```

![](Diversity_Analysis_files/figure-gfm/unnamed-chunk-28-1.png)<!-- -->

``` r
dist_site_tukey <- TukeyHSD(metagbeta_clr_site)
microbe_site_tuk <- data.frame(dist_site_tukey[['group']])
write.csv(microbe_site_tuk, "./clr_betadispr_site_tukey.csv")
ggsave('./clr_betadispr_box_site.png',width=12, height=4, plot = dist_microbe_box_site)
```

``` r
dist_microbe_box_date <- ggboxplot(metagdistances_date, x = "date", y = "Distance_to_Centroid" , add = "point")
dist_microbe_box_date
```

![](Diversity_Analysis_files/figure-gfm/unnamed-chunk-29-1.png)<!-- -->

``` r
dist_date_tukey <- TukeyHSD(metagbeta_clr_date)
microbe_date_tuk <- data.frame(dist_date_tukey[['group']])
write.csv(microbe_date_tuk, "./clr_betadispr_date_tukey.csv")
ggsave('./clr_betadispr_box_date.png',width=12, height=4, plot = dist_microbe_box_date)
```

``` r
betadipsr_figure_site_date <- ggarrange(dist_microbe_box_site,dist_microbe_box_date,
                                    labels = c("A", "B"),
                                    ncol = 1, nrow = 2)
betadipsr_figure_site_date
```

![](Diversity_Analysis_files/figure-gfm/unnamed-chunk-30-1.png)<!-- -->

``` r
ggsave('./clr_betadispr_boxplot_site_date.png',width=15, height=8, plot = betadipsr_figure_site_date)
```
