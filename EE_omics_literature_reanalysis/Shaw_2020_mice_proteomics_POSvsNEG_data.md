2025_10_28_Shaw_2020_proteomics_data
================
Stephanie Huang
2025-12-30

NOTE:

Exported file only includes Rats & mouse mapping (excluded the 210 other
rodents)

Also, the provided data set only seems to have A1-A3 (in the frontal) &
B1-B3 (in hippo file). Genuinely super confused by their data… from the
looks of it it seem like like they did ~5 technical replicates per
sample.

Links

- Dataset via Mendeley Data -
  <https://data.mendeley.com/datasets/5b4477386w/1>

Legit kinda weird! - Has a billion rodents, the mapping of groups &
samples is confusing

# Part 1 - Setup & General Checks

Part 1 contains the setup & general checks. This data is really weird so
here I make various plots & tables to examine the data.

## Setup

``` r
library(tidyverse) # ggplot, dplyr, etc
library(ggrepel) #  for geom_text_repel() on volcanos
library(gt) # gt() tables
library(pheatmap) # heatmaps
library(ggh4x) # facet_nested
library(readxl)
library(rmarkdown) # already loaded or no?? -> for github markdown test
```

### Load in data

``` r
frontal_file <- read_csv("input/Shaw_2020_UCLan_Preomics_Data_Frontal.csv", col_names = FALSE)
```

    ## Rows: 980 Columns: 28
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: ","
    ## chr (28): X1, X2, X3, X4, X5, X6, X7, X8, X9, X10, X11, X12, X13, X14, X15, ...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
colnames(frontal_file) <- frontal_file[3,] # assign row 3 as headers

frontal_file <- frontal_file[-1:-3,] # negative index (rm rows 1-3) |data starts r4

# 1st data obs is ODP2_MESAU 
```

``` r
head(frontal_file)
```

    ## # A tibble: 6 × 28
    ##   Accession     `Peptide count` `Unique peptides` `Confidence score` `Anova (p)`
    ##   <chr>         <chr>           <chr>             <chr>              <chr>      
    ## 1 ODP2_MESAU    7               0                 161.8              0.002464258
    ## 2 STAG3_RAT     4               1                 15.18              0.00260158 
    ## 3 ZFP11_MOUSE   3               2                 33.49              0.003140858
    ## 4 KPCA_MOUSE;K… 10              3                 131.2              0.003190691
    ## 5 APT_MOUSE     2               0                 41.77              0.003780238
    ## 6 DYN3_RAT;DYN… 11              1                 259.53             0.006541435
    ## # ℹ 23 more variables: `q Value` <chr>, `Max fold change` <chr>, Power <chr>,
    ## #   `Highest mean condition` <chr>, `Lowest mean condition` <chr>, Mass <chr>,
    ## #   Description <chr>, AS_01_1_in_10 <chr>, AS_01_1_in_10_5ul <chr>,
    ## #   AS_1_in_10_03 <chr>, AS_1_in_10_25 <chr>, AS_1_in_10_27 <chr>,
    ## #   AS_1_in_10_29 <chr>, AS_1_in_10_05 <chr>, AS_1_in_10_07 <chr>,
    ## #   AS_1_in_10_09 <chr>, AS_1_in_10_11 <chr>, AS_1_in_10_13 <chr>,
    ## #   AS_1_in_10_15 <chr>, AS_1_in_10_17 <chr>, AS_1_in_10_19 <chr>, …

``` r
str(frontal_file)
```

    ## tibble [977 × 28] (S3: tbl_df/tbl/data.frame)
    ##  $ Accession             : chr [1:977] "ODP2_MESAU" "STAG3_RAT" "ZFP11_MOUSE" "KPCA_MOUSE;KPCA_RAT" ...
    ##  $ Peptide count         : chr [1:977] "7" "4" "3" "10" ...
    ##  $ Unique peptides       : chr [1:977] "0" "1" "2" "3" ...
    ##  $ Confidence score      : chr [1:977] "161.8" "15.18" "33.49" "131.2" ...
    ##  $ Anova (p)             : chr [1:977] "0.002464258" "0.00260158" "0.003140858" "0.003190691" ...
    ##  $ q Value               : chr [1:977] "0.023692053" "0.023692053" "0.023692053" "0.023692053" ...
    ##  $ Max fold change       : chr [1:977] "1.601685259" "1.952788434" "2.52415046" "3.783920472" ...
    ##  $ Power                 : chr [1:977] "0.944748878" "0.935873573" "0.938307423" "0.933402809" ...
    ##  $ Highest mean condition: chr [1:977] "2A" "3A" "2A" "2A" ...
    ##  $ Lowest mean condition : chr [1:977] "1A" "1A" "1A" "1A" ...
    ##  $ Mass                  : chr [1:977] "23273" "143628" "78691" "77943" ...
    ##  $ Description           : chr [1:977] "Dihydrolipoyllysine-residue acetyltransferase component of pyruvate dehydrogenase complex, mitochondrial (Fragm"| __truncated__ "Cohesin subunit SA-3 OS=Rattus norvegicus GN=Stag3 PE=2 SV=1" "Zinc finger protein 11 OS=Mus musculus GN=Zfp11 PE=2 SV=2" "Protein kinase C alpha type OS=Mus musculus GN=Prkca PE=1 SV=3" ...
    ##  $ AS_01_1_in_10         : chr [1:977] "1536900.418" "713004.757" "343529.7736" "1111827.605" ...
    ##  $ AS_01_1_in_10_5ul     : chr [1:977] "2569006.718" "1645355.064" "479218.4743" "3182320.255" ...
    ##  $ AS_1_in_10_03         : chr [1:977] "2448567.374" "1560395.28" "348831.7091" "2727825.211" ...
    ##  $ AS_1_in_10_25         : chr [1:977] "2069393.943" "1155522.648" "462370.9874" "4959118.491" ...
    ##  $ AS_1_in_10_27         : chr [1:977] "1969278.399" "1444907.75" "443755.6015" "1657414.183" ...
    ##  $ AS_1_in_10_29         : chr [1:977] "2650719.779" "1519601.789" "311857.359" "2285505.521" ...
    ##  $ AS_1_in_10_05         : chr [1:977] "3404033.334" "1842621.712" "804864.0543" "5970597.607" ...
    ##  $ AS_1_in_10_07         : chr [1:977] "3936013.466" "3720515.192" "1077451.648" "15392801.54" ...
    ##  $ AS_1_in_10_09         : chr [1:977] "4085866.87" "2294235.372" "1993207.108" "18083033.44" ...
    ##  $ AS_1_in_10_11         : chr [1:977] "2662270.871" "2328160.147" "490472.1177" "4006316.076" ...
    ##  $ AS_1_in_10_13         : chr [1:977] "3588903.754" "2418806.469" "660354.0951" "6759911.527" ...
    ##  $ AS_1_in_10_15         : chr [1:977] "2642625.201" "2749643.298" "486490.8217" "4053437.83" ...
    ##  $ AS_1_in_10_17         : chr [1:977] "3080639.923" "2796708.291" "566600.2184" "4866733.059" ...
    ##  $ AS_1_in_10_19         : chr [1:977] "2238860.372" "3709033.187" "433968.5678" "3030507.848" ...
    ##  $ AS_1_in_10_21         : chr [1:977] "3232330.752" "1961915.98" "450165.0993" "4135319.668" ...
    ##  $ AS_1_in_10_23         : chr [1:977] "2769823.436" "1864408.278" "478286.8019" "3537549.281" ...

``` r
hippo_file <- read_csv("input/Shaw_2020_UCLan_Preomics_Data_Hippocampus.csv", col_names = FALSE)
```

    ## Rows: 980 Columns: 27
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: ","
    ## chr (27): X1, X2, X3, X4, X5, X6, X7, X8, X9, X10, X11, X12, X13, X14, X15, ...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
colnames(hippo_file) <- hippo_file[3,] # assign row 3 as header

hippo_file <- hippo_file[-1:-3,] # negative index (rm rows 1-3) |data starts r4
# 1st data obs is KPCA mouse
```

``` r
head(hippo_file)
```

    ## # A tibble: 6 × 27
    ##   Accession     `Peptide count` `Unique peptides` `Confidence score` `Anova (p)`
    ##   <chr>         <chr>           <chr>             <chr>              <chr>      
    ## 1 KPCA_MOUSE;K… 10              3                 131.2              0.002910233
    ## 2 UBE2N_MOUSE   3               1                 110.32             0.017927874
    ## 3 CPNE3_MOUSE   2               0                 5.24               0.025786986
    ## 4 CMC1_MOUSE    10              3                 383.05             0.032557404
    ## 5 ATPD_MOUSE    5               2                 173.48             0.046283841
    ## 6 SUMO3_MOUSE   2               0                 48.17              0.051328313
    ## # ℹ 22 more variables: `q Value` <chr>, `Max fold change` <chr>, Power <chr>,
    ## #   `Highest mean condition` <chr>, `Lowest mean condition` <chr>, Mass <chr>,
    ## #   Description <chr>, AS_1_in_10_02 <chr>, AS_1_in_10_04 <chr>,
    ## #   AS_1_in_10_26 <chr>, AS_1_in_10_28 <chr>, AS_1_in_10_30 <chr>,
    ## #   AS_1_in_10_06 <chr>, AS_1_in_10_08 <chr>, AS_1_in_10_10 <chr>,
    ## #   AS_1_in_10_12 <chr>, AS_1_in_10_14 <chr>, AS_1_in_10_16 <chr>,
    ## #   AS_1_in_10_18 <chr>, AS_1_in_10_20 <chr>, AS_1_in_10_22 <chr>, …

``` r
str(hippo_file)
```

    ## tibble [977 × 27] (S3: tbl_df/tbl/data.frame)
    ##  $ Accession             : chr [1:977] "KPCA_MOUSE;KPCA_RAT" "UBE2N_MOUSE" "CPNE3_MOUSE" "CMC1_MOUSE" ...
    ##  $ Peptide count         : chr [1:977] "10" "3" "2" "10" ...
    ##  $ Unique peptides       : chr [1:977] "3" "1" "0" "3" ...
    ##  $ Confidence score      : chr [1:977] "131.2" "110.32" "5.24" "383.05" ...
    ##  $ Anova (p)             : chr [1:977] "0.002910233" "0.017927874" "0.025786986" "0.032557404" ...
    ##  $ q Value               : chr [1:977] "0.999952175" "0.999952175" "0.999952175" "0.999952175" ...
    ##  $ Max fold change       : chr [1:977] "3.683703516" "2.935722795" "2.709096029" "1.60390951" ...
    ##  $ Power                 : chr [1:977] "0.945790397" "0.762659981" "0.705040749" "0.664884067" ...
    ##  $ Highest mean condition: chr [1:977] "2B" "2B" "1B" "3B" ...
    ##  $ Lowest mean condition : chr [1:977] "1B" "1B" "3B" "1B" ...
    ##  $ Mass                  : chr [1:977] "77943" "17184" "60288" "74922" ...
    ##  $ Description           : chr [1:977] "Protein kinase C alpha type OS=Mus musculus GN=Prkca PE=1 SV=3" "Ubiquitin-conjugating enzyme E2 N OS=Mus musculus GN=Ube2n PE=1 SV=1" "Copine-3 OS=Mus musculus GN=Cpne3 PE=1 SV=2" "Calcium-binding mitochondrial carrier protein Aralar1 OS=Mus musculus GN=Slc25a12 PE=1 SV=1" ...
    ##  $ AS_1_in_10_02         : chr [1:977] "2951760.5" "157193.3077" "947689.8367" "588198.9334" ...
    ##  $ AS_1_in_10_04         : chr [1:977] "2291970.901" "91989.11491" "1319530.082" "332668.8275" ...
    ##  $ AS_1_in_10_26         : chr [1:977] "923192.5048" "206082.7776" "1245403.509" "591568.7441" ...
    ##  $ AS_1_in_10_28         : chr [1:977] "3135994.946" "163809.9169" "603469.4988" "581131.6073" ...
    ##  $ AS_1_in_10_30         : chr [1:977] "2437314.871" "129828.4942" "511829.6218" "523388.5588" ...
    ##  $ AS_1_in_10_06         : chr [1:977] "5872546.785" "233327.4842" "797742.6897" "678359.6556" ...
    ##  $ AS_1_in_10_08         : chr [1:977] "7806426.002" "667819.9435" "1107893.334" "949742.7586" ...
    ##  $ AS_1_in_10_10         : chr [1:977] "16421694.6" "756964.1161" "260543.2857" "672135.2727" ...
    ##  $ AS_1_in_10_12         : chr [1:977] "7096958.665" "162271.0907" "642863.4316" "893978.7535" ...
    ##  $ AS_1_in_10_14         : chr [1:977] "6049914.193" "378190.7684" "326191.4471" "783372.7912" ...
    ##  $ AS_1_in_10_16         : chr [1:977] "6369199.382" "244331.0006" "321523.6919" "863131.1426" ...
    ##  $ AS_1_in_10_18         : chr [1:977] "12260291.73" "204360.127" "231196.5845" "1108353.077" ...
    ##  $ AS_1_in_10_20         : chr [1:977] "5850381.165" "121085.4446" "178237.6543" "1108832.39" ...
    ##  $ AS_1_in_10_22         : chr [1:977] "3904894.457" "98114.30979" "282212.0763" "514850.6123" ...
    ##  $ AS_1_in_10_24         : chr [1:977] "3176707.808" "212338.5481" "695120.3223" "602194.4697" ...

``` r
fr_others <- frontal_file %>%
  filter(!str_detect(Accession, "MOUSE|RAT")) # not mouse or rat

# gives weird rodents - 80 of 977 are not mice or rats 
# like GNAI2_CAVPO (Gerbil) ODP2_MESAU Golden Hamster | GNAI2_CAVPO Guinea Pig
```

### Clean up & bind

(doing it sep as the col headers are slightly different)

``` r
# frontal headers
frontal_headers <- data.frame(
  col_ids = c('AS_01_1_in_10',  'AS_01_1_in_10_5ul',    'AS_1_in_10_03',
                 'AS_1_in_10_25',   'AS_1_in_10_27',    'AS_1_in_10_29',
                 'AS_1_in_10_05',   'AS_1_in_10_07',    'AS_1_in_10_09',
                 'AS_1_in_10_11',   'AS_1_in_10_13',    'AS_1_in_10_15',
                 'AS_1_in_10_17',   'AS_1_in_10_19',    'AS_1_in_10_21',
                 'AS_1_in_10_23'),
  col_headers = c(rep('1A', 6), rep('2A', 5), rep('3A', 5))
  )

# frontal headers
hip_headers <- data.frame(
  col_ids = c('AS_1_in_10_02',  'AS_1_in_10_04',    'AS_1_in_10_26',
              'AS_1_in_10_28',  'AS_1_in_10_30',    'AS_1_in_10_06',
              'AS_1_in_10_08',  'AS_1_in_10_10',    'AS_1_in_10_12',
              'AS_1_in_10_14',  'AS_1_in_10_16',    'AS_1_in_10_18',
              'AS_1_in_10_20',  'AS_1_in_10_22',    'AS_1_in_10_24'),
  col_headers = c(rep('1B', 5), rep('2B', 5), rep('3B', 5)))
```

``` r
# fr_mice <- frontal_file %>%
#   filter(str_detect(Accession, "MOUSE"))
# 
# # 675/977 contain mice, 
```

``` r
Frontal <- frontal_file %>%
  # filter(`Anova (p)` < 0.0505) %>%
  select(1:12) %>%
  mutate_at(c('Confidence score', 'Anova (p)', 'q Value',
              'Max fold change', 'Power'), as.numeric)
```

``` r
Hippo <- hippo_file %>%
  # filter(`Anova (p)` < 0.0505) %>%
  select(1:12) %>%
  mutate_at(c('Confidence score', 'Anova (p)', 'q Value',
              'Max fold change', 'Power'), as.numeric)
```

``` r
combined <- bind_rows(list(Frontal = Frontal, Hippocampus = Hippo),
                      .id = "Region") %>%
  mutate( # recode so 1A 2B match the groups - A is POS, B is NEG
    highest_mean = ifelse(str_detect(`Highest mean condition`, "A"), "POS", "NEG"),
    lowest_mean = ifelse(str_detect(`Lowest mean condition`, "A"), "POS", "NEG"),
    .after = `Lowest mean condition`
    )
```

``` r
summary(combined)
```

    ##     Region           Accession         Peptide count      Unique peptides   
    ##  Length:1954        Length:1954        Length:1954        Length:1954       
    ##  Class :character   Class :character   Class :character   Class :character  
    ##  Mode  :character   Mode  :character   Mode  :character   Mode  :character  
    ##                                                                             
    ##                                                                             
    ##                                                                             
    ##  Confidence score    Anova (p)           q Value        Max fold change  
    ##  Min.   :   1.28   Min.   :0.002464   Min.   :0.02369   Min.   :  1.005  
    ##  1st Qu.:  61.55   1st Qu.:0.149345   1st Qu.:0.02369   1st Qu.:  1.353  
    ##  Median : 114.50   Median :0.326242   Median :0.54031   Median :  1.714  
    ##  Mean   : 237.07   Mean   :0.446985   Mean   :0.51327   Mean   :  1.931  
    ##  3rd Qu.: 232.81   3rd Qu.:0.776132   3rd Qu.:0.99995   3rd Qu.:  2.110  
    ##  Max.   :3215.49   Max.   :0.999952   Max.   :0.99995   Max.   :123.442  
    ##      Power         Highest mean condition Lowest mean condition
    ##  Min.   :0.05001   Length:1954            Length:1954          
    ##  1st Qu.:0.08190   Class :character       Class :character     
    ##  Median :0.21717   Mode  :character       Mode  :character     
    ##  Mean   :0.23555                                               
    ##  3rd Qu.:0.35831                                               
    ##  Max.   :0.94579                                               
    ##  highest_mean       lowest_mean            Mass           Description       
    ##  Length:1954        Length:1954        Length:1954        Length:1954       
    ##  Class :character   Class :character   Class :character   Class :character  
    ##  Mode  :character   Mode  :character   Mode  :character   Mode  :character  
    ##                                                                             
    ##                                                                             
    ## 

``` r
volcano_pv <- ggplot(combined, 
                     aes(x = `Max fold change`, y = -log10(`Anova (p)`),
                         label = Accession)) +
  geom_point(aes(color = case_when(
    `Anova (p)` < 0.05 & `Max fold change` > 0.05 ~ "Upregulated (p<0.05)",
    `Anova (p)` < 0.05 & `Max fold change` < -0.05 ~ "Downregulated (p<0.05)",
                                   TRUE ~ "Not Significant")),
            alpha = 1, size = 2, shape = 1) +
  xlim(-8,12) +
  geom_vline(xintercept = 0, linetype = "dashed",
             color = "grey40", linewidth = .5) +
  scale_color_manual(values = c("firebrick3", 
                                "royalblue4", 
                                 "grey80"),
                     breaks = c("Upregulated (p<0.05)", 
                                "Downregulated (p<0.05)", 
                                "Not Significant"),
                     labels = c("Upregulated Proteins (p<0.05)", 
                                "Downregulated Proteins (p<0.05)", 
                                "Not Significant")) +
  
  geom_text_repel() +
  facet_wrap(highest_mean~Region, ncol = 2) +
  labs(
    title = "???? no fold change? or is this just up in EE?",
    subtitle = 'Also coding is weird - why frontal all higes in pos & neg',
       x = "Log2 Signal",
       y = "-Log10(p-value)",
       color = "mRNA Regulation") +
  theme_bw() +
  theme(legend.position = "bottom")



volcano_pv
```

    ## Warning: Removed 4 rows containing missing values or values outside the scale range
    ## (`geom_point()`).

    ## Warning: Removed 4 rows containing missing values or values outside the scale range
    ## (`geom_text_repel()`).

    ## Warning: ggrepel: 965 unlabeled data points (too many overlaps). Consider
    ## increasing max.overlaps

    ## Warning: ggrepel: 966 unlabeled data points (too many overlaps). Consider
    ## increasing max.overlaps

![](Shaw_2020_mice_proteomics_POSvsNEG_data_files/figure-gfm/unnamed-chunk-13-1.png)<!-- -->

------------------------------------------------------------------------

### Heatmap check (done on separate raw/original df provided by shaw)

``` r
# Making a maxtrix - dplyr/tidyverse style
f_matrix <- frontal_file %>%
  filter(`Anova (p)` <= 0.06) %>%
  select(c(Accession,AS_01_1_in_10:AS_1_in_10_23)) %>% # get cols
  tibble::column_to_rownames("Accession") %>% # set Accession to rowname 
  mutate(across(everything(), as.numeric)) %>% # as numeric! 
  as.matrix() # make matrix 

# 
f_heatmap <- pheatmap(f_matrix, main = "Frontal Lobe - p.05")
```

![](Shaw_2020_mice_proteomics_POSvsNEG_data_files/figure-gfm/unnamed-chunk-14-1.png)<!-- -->

``` r
# Making a maxtrix - dplyr/tidyverse style
h_matrix <- hippo_file %>%
  filter(`Anova (p)` <= 0.05) %>%
  select(c(Accession,AS_1_in_10_02:AS_1_in_10_24)) %>% # get cols
  tibble::column_to_rownames("Accession") %>% # set Accession to rowname
  mutate(across(everything(), as.numeric)) %>% # as numeric!
  as.matrix() # make matrix

h_heatmap <- pheatmap(h_matrix, main = "Hippocampus - p.05")
```

![](Shaw_2020_mice_proteomics_POSvsNEG_data_files/figure-gfm/unnamed-chunk-15-1.png)<!-- -->

------------------------------------------------------------------------

### gt table (of sig)

``` r
combined |>
  dplyr::filter(`Anova (p)` <= 0.05) |>
  gt(
    groupname_col = "Region"
    ) |>
  fmt_number(decimals = 3) |>
  # data_color(columns = "Regulation") +
  tab_options(table.font.size=11.5) |>
  tab_header(
    title = "Shaw et al., 2020 Proteomics - NEG vs POS - Check 1",
    subtitle = "From Frontal (47/48) & Hippocampus (5) - p-value <0.0505"
  )
```

<div id="tkcaxwsnyc" style="padding-left:0px;padding-right:0px;padding-top:10px;padding-bottom:10px;overflow-x:auto;overflow-y:auto;width:auto;height:auto;">
<style>#tkcaxwsnyc table {
  font-family: system-ui, 'Segoe UI', Roboto, Helvetica, Arial, sans-serif, 'Apple Color Emoji', 'Segoe UI Emoji', 'Segoe UI Symbol', 'Noto Color Emoji';
  -webkit-font-smoothing: antialiased;
  -moz-osx-font-smoothing: grayscale;
}
&#10;#tkcaxwsnyc thead, #tkcaxwsnyc tbody, #tkcaxwsnyc tfoot, #tkcaxwsnyc tr, #tkcaxwsnyc td, #tkcaxwsnyc th {
  border-style: none;
}
&#10;#tkcaxwsnyc p {
  margin: 0;
  padding: 0;
}
&#10;#tkcaxwsnyc .gt_table {
  display: table;
  border-collapse: collapse;
  line-height: normal;
  margin-left: auto;
  margin-right: auto;
  color: #333333;
  font-size: 11.5px;
  font-weight: normal;
  font-style: normal;
  background-color: #FFFFFF;
  width: auto;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #A8A8A8;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #A8A8A8;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
}
&#10;#tkcaxwsnyc .gt_caption {
  padding-top: 4px;
  padding-bottom: 4px;
}
&#10;#tkcaxwsnyc .gt_title {
  color: #333333;
  font-size: 125%;
  font-weight: initial;
  padding-top: 4px;
  padding-bottom: 4px;
  padding-left: 5px;
  padding-right: 5px;
  border-bottom-color: #FFFFFF;
  border-bottom-width: 0;
}
&#10;#tkcaxwsnyc .gt_subtitle {
  color: #333333;
  font-size: 85%;
  font-weight: initial;
  padding-top: 3px;
  padding-bottom: 5px;
  padding-left: 5px;
  padding-right: 5px;
  border-top-color: #FFFFFF;
  border-top-width: 0;
}
&#10;#tkcaxwsnyc .gt_heading {
  background-color: #FFFFFF;
  text-align: center;
  border-bottom-color: #FFFFFF;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
}
&#10;#tkcaxwsnyc .gt_bottom_border {
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}
&#10;#tkcaxwsnyc .gt_col_headings {
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
}
&#10;#tkcaxwsnyc .gt_col_heading {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: normal;
  text-transform: inherit;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: bottom;
  padding-top: 5px;
  padding-bottom: 6px;
  padding-left: 5px;
  padding-right: 5px;
  overflow-x: hidden;
}
&#10;#tkcaxwsnyc .gt_column_spanner_outer {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: normal;
  text-transform: inherit;
  padding-top: 0;
  padding-bottom: 0;
  padding-left: 4px;
  padding-right: 4px;
}
&#10;#tkcaxwsnyc .gt_column_spanner_outer:first-child {
  padding-left: 0;
}
&#10;#tkcaxwsnyc .gt_column_spanner_outer:last-child {
  padding-right: 0;
}
&#10;#tkcaxwsnyc .gt_column_spanner {
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  vertical-align: bottom;
  padding-top: 5px;
  padding-bottom: 5px;
  overflow-x: hidden;
  display: inline-block;
  width: 100%;
}
&#10;#tkcaxwsnyc .gt_spanner_row {
  border-bottom-style: hidden;
}
&#10;#tkcaxwsnyc .gt_group_heading {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  text-transform: inherit;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: middle;
  text-align: left;
}
&#10;#tkcaxwsnyc .gt_empty_group_heading {
  padding: 0.5px;
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  vertical-align: middle;
}
&#10;#tkcaxwsnyc .gt_from_md > :first-child {
  margin-top: 0;
}
&#10;#tkcaxwsnyc .gt_from_md > :last-child {
  margin-bottom: 0;
}
&#10;#tkcaxwsnyc .gt_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  margin: 10px;
  border-top-style: solid;
  border-top-width: 1px;
  border-top-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: middle;
  overflow-x: hidden;
}
&#10;#tkcaxwsnyc .gt_stub {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  text-transform: inherit;
  border-right-style: solid;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
  padding-left: 5px;
  padding-right: 5px;
}
&#10;#tkcaxwsnyc .gt_stub_row_group {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  text-transform: inherit;
  border-right-style: solid;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
  padding-left: 5px;
  padding-right: 5px;
  vertical-align: top;
}
&#10;#tkcaxwsnyc .gt_row_group_first td {
  border-top-width: 2px;
}
&#10;#tkcaxwsnyc .gt_row_group_first th {
  border-top-width: 2px;
}
&#10;#tkcaxwsnyc .gt_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}
&#10;#tkcaxwsnyc .gt_first_summary_row {
  border-top-style: solid;
  border-top-color: #D3D3D3;
}
&#10;#tkcaxwsnyc .gt_first_summary_row.thick {
  border-top-width: 2px;
}
&#10;#tkcaxwsnyc .gt_last_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}
&#10;#tkcaxwsnyc .gt_grand_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}
&#10;#tkcaxwsnyc .gt_first_grand_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-top-style: double;
  border-top-width: 6px;
  border-top-color: #D3D3D3;
}
&#10;#tkcaxwsnyc .gt_last_grand_summary_row_top {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-bottom-style: double;
  border-bottom-width: 6px;
  border-bottom-color: #D3D3D3;
}
&#10;#tkcaxwsnyc .gt_striped {
  background-color: rgba(128, 128, 128, 0.05);
}
&#10;#tkcaxwsnyc .gt_table_body {
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}
&#10;#tkcaxwsnyc .gt_footnotes {
  color: #333333;
  background-color: #FFFFFF;
  border-bottom-style: none;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
}
&#10;#tkcaxwsnyc .gt_footnote {
  margin: 0px;
  font-size: 90%;
  padding-top: 4px;
  padding-bottom: 4px;
  padding-left: 5px;
  padding-right: 5px;
}
&#10;#tkcaxwsnyc .gt_sourcenotes {
  color: #333333;
  background-color: #FFFFFF;
  border-bottom-style: none;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
}
&#10;#tkcaxwsnyc .gt_sourcenote {
  font-size: 90%;
  padding-top: 4px;
  padding-bottom: 4px;
  padding-left: 5px;
  padding-right: 5px;
}
&#10;#tkcaxwsnyc .gt_left {
  text-align: left;
}
&#10;#tkcaxwsnyc .gt_center {
  text-align: center;
}
&#10;#tkcaxwsnyc .gt_right {
  text-align: right;
  font-variant-numeric: tabular-nums;
}
&#10;#tkcaxwsnyc .gt_font_normal {
  font-weight: normal;
}
&#10;#tkcaxwsnyc .gt_font_bold {
  font-weight: bold;
}
&#10;#tkcaxwsnyc .gt_font_italic {
  font-style: italic;
}
&#10;#tkcaxwsnyc .gt_super {
  font-size: 65%;
}
&#10;#tkcaxwsnyc .gt_footnote_marks {
  font-size: 75%;
  vertical-align: 0.4em;
  position: initial;
}
&#10;#tkcaxwsnyc .gt_asterisk {
  font-size: 100%;
  vertical-align: 0;
}
&#10;#tkcaxwsnyc .gt_indent_1 {
  text-indent: 5px;
}
&#10;#tkcaxwsnyc .gt_indent_2 {
  text-indent: 10px;
}
&#10;#tkcaxwsnyc .gt_indent_3 {
  text-indent: 15px;
}
&#10;#tkcaxwsnyc .gt_indent_4 {
  text-indent: 20px;
}
&#10;#tkcaxwsnyc .gt_indent_5 {
  text-indent: 25px;
}
&#10;#tkcaxwsnyc .katex-display {
  display: inline-flex !important;
  margin-bottom: 0.75em !important;
}
&#10;#tkcaxwsnyc div.Reactable > div.rt-table > div.rt-thead > div.rt-tr.rt-tr-group-header > div.rt-th-group:after {
  height: 0px !important;
}
</style>
<table class="gt_table" data-quarto-disable-processing="false" data-quarto-bootstrap="false">
  <thead>
    <tr class="gt_heading">
      <td colspan="14" class="gt_heading gt_title gt_font_normal" style>Shaw et al., 2020 Proteomics - NEG vs POS - Check 1</td>
    </tr>
    <tr class="gt_heading">
      <td colspan="14" class="gt_heading gt_subtitle gt_font_normal gt_bottom_border" style>From Frontal (47/48) &amp; Hippocampus (5) - p-value &lt;0.0505</td>
    </tr>
    <tr class="gt_col_headings">
      <th class="gt_col_heading gt_columns_bottom_border gt_left" rowspan="1" colspan="1" scope="col" id="Accession">Accession</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="Peptide-count">Peptide count</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="Unique-peptides">Unique peptides</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="Confidence-score">Confidence score</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="Anova-(p)">Anova (p)</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="q-Value">q Value</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="Max-fold-change">Max fold change</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="Power">Power</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_left" rowspan="1" colspan="1" scope="col" id="Highest-mean-condition">Highest mean condition</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_left" rowspan="1" colspan="1" scope="col" id="Lowest-mean-condition">Lowest mean condition</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_left" rowspan="1" colspan="1" scope="col" id="highest_mean">highest_mean</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_left" rowspan="1" colspan="1" scope="col" id="lowest_mean">lowest_mean</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="Mass">Mass</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_left" rowspan="1" colspan="1" scope="col" id="Description">Description</th>
    </tr>
  </thead>
  <tbody class="gt_table_body">
    <tr class="gt_group_heading_row">
      <th colspan="14" class="gt_group_heading" scope="colgroup" id="Frontal">Frontal</th>
    </tr>
    <tr class="gt_row_group_first"><td headers="Frontal  Accession" class="gt_row gt_left">ODP2_MESAU</td>
<td headers="Frontal  Peptide count" class="gt_row gt_right">7</td>
<td headers="Frontal  Unique peptides" class="gt_row gt_right">0</td>
<td headers="Frontal  Confidence score" class="gt_row gt_right">161.800</td>
<td headers="Frontal  Anova (p)" class="gt_row gt_right">0.002</td>
<td headers="Frontal  q Value" class="gt_row gt_right">0.024</td>
<td headers="Frontal  Max fold change" class="gt_row gt_right">1.602</td>
<td headers="Frontal  Power" class="gt_row gt_right">0.945</td>
<td headers="Frontal  Highest mean condition" class="gt_row gt_left">2A</td>
<td headers="Frontal  Lowest mean condition" class="gt_row gt_left">1A</td>
<td headers="Frontal  highest_mean" class="gt_row gt_left">POS</td>
<td headers="Frontal  lowest_mean" class="gt_row gt_left">POS</td>
<td headers="Frontal  Mass" class="gt_row gt_right">23273</td>
<td headers="Frontal  Description" class="gt_row gt_left">Dihydrolipoyllysine-residue acetyltransferase component of pyruvate dehydrogenase complex, mitochondrial (Fragments) OS=Mesocricetus auratus GN=DLAT PE=1 SV=1</td></tr>
    <tr><td headers="Frontal  Accession" class="gt_row gt_left">STAG3_RAT</td>
<td headers="Frontal  Peptide count" class="gt_row gt_right">4</td>
<td headers="Frontal  Unique peptides" class="gt_row gt_right">1</td>
<td headers="Frontal  Confidence score" class="gt_row gt_right">15.180</td>
<td headers="Frontal  Anova (p)" class="gt_row gt_right">0.003</td>
<td headers="Frontal  q Value" class="gt_row gt_right">0.024</td>
<td headers="Frontal  Max fold change" class="gt_row gt_right">1.953</td>
<td headers="Frontal  Power" class="gt_row gt_right">0.936</td>
<td headers="Frontal  Highest mean condition" class="gt_row gt_left">3A</td>
<td headers="Frontal  Lowest mean condition" class="gt_row gt_left">1A</td>
<td headers="Frontal  highest_mean" class="gt_row gt_left">POS</td>
<td headers="Frontal  lowest_mean" class="gt_row gt_left">POS</td>
<td headers="Frontal  Mass" class="gt_row gt_right">143628</td>
<td headers="Frontal  Description" class="gt_row gt_left">Cohesin subunit SA-3 OS=Rattus norvegicus GN=Stag3 PE=2 SV=1</td></tr>
    <tr><td headers="Frontal  Accession" class="gt_row gt_left">ZFP11_MOUSE</td>
<td headers="Frontal  Peptide count" class="gt_row gt_right">3</td>
<td headers="Frontal  Unique peptides" class="gt_row gt_right">2</td>
<td headers="Frontal  Confidence score" class="gt_row gt_right">33.490</td>
<td headers="Frontal  Anova (p)" class="gt_row gt_right">0.003</td>
<td headers="Frontal  q Value" class="gt_row gt_right">0.024</td>
<td headers="Frontal  Max fold change" class="gt_row gt_right">2.524</td>
<td headers="Frontal  Power" class="gt_row gt_right">0.938</td>
<td headers="Frontal  Highest mean condition" class="gt_row gt_left">2A</td>
<td headers="Frontal  Lowest mean condition" class="gt_row gt_left">1A</td>
<td headers="Frontal  highest_mean" class="gt_row gt_left">POS</td>
<td headers="Frontal  lowest_mean" class="gt_row gt_left">POS</td>
<td headers="Frontal  Mass" class="gt_row gt_right">78691</td>
<td headers="Frontal  Description" class="gt_row gt_left">Zinc finger protein 11 OS=Mus musculus GN=Zfp11 PE=2 SV=2</td></tr>
    <tr><td headers="Frontal  Accession" class="gt_row gt_left">KPCA_MOUSE;KPCA_RAT</td>
<td headers="Frontal  Peptide count" class="gt_row gt_right">10</td>
<td headers="Frontal  Unique peptides" class="gt_row gt_right">3</td>
<td headers="Frontal  Confidence score" class="gt_row gt_right">131.200</td>
<td headers="Frontal  Anova (p)" class="gt_row gt_right">0.003</td>
<td headers="Frontal  q Value" class="gt_row gt_right">0.024</td>
<td headers="Frontal  Max fold change" class="gt_row gt_right">3.784</td>
<td headers="Frontal  Power" class="gt_row gt_right">0.933</td>
<td headers="Frontal  Highest mean condition" class="gt_row gt_left">2A</td>
<td headers="Frontal  Lowest mean condition" class="gt_row gt_left">1A</td>
<td headers="Frontal  highest_mean" class="gt_row gt_left">POS</td>
<td headers="Frontal  lowest_mean" class="gt_row gt_left">POS</td>
<td headers="Frontal  Mass" class="gt_row gt_right">77943</td>
<td headers="Frontal  Description" class="gt_row gt_left">Protein kinase C alpha type OS=Mus musculus GN=Prkca PE=1 SV=3</td></tr>
    <tr><td headers="Frontal  Accession" class="gt_row gt_left">APT_MOUSE</td>
<td headers="Frontal  Peptide count" class="gt_row gt_right">2</td>
<td headers="Frontal  Unique peptides" class="gt_row gt_right">0</td>
<td headers="Frontal  Confidence score" class="gt_row gt_right">41.770</td>
<td headers="Frontal  Anova (p)" class="gt_row gt_right">0.004</td>
<td headers="Frontal  q Value" class="gt_row gt_right">0.024</td>
<td headers="Frontal  Max fold change" class="gt_row gt_right">3.507</td>
<td headers="Frontal  Power" class="gt_row gt_right">0.917</td>
<td headers="Frontal  Highest mean condition" class="gt_row gt_left">3A</td>
<td headers="Frontal  Lowest mean condition" class="gt_row gt_left">1A</td>
<td headers="Frontal  highest_mean" class="gt_row gt_left">POS</td>
<td headers="Frontal  lowest_mean" class="gt_row gt_left">POS</td>
<td headers="Frontal  Mass" class="gt_row gt_right">19883</td>
<td headers="Frontal  Description" class="gt_row gt_left">Adenine phosphoribosyltransferase OS=Mus musculus GN=Aprt PE=1 SV=2</td></tr>
    <tr><td headers="Frontal  Accession" class="gt_row gt_left">DYN3_RAT;DYN3_MOUSE</td>
<td headers="Frontal  Peptide count" class="gt_row gt_right">11</td>
<td headers="Frontal  Unique peptides" class="gt_row gt_right">1</td>
<td headers="Frontal  Confidence score" class="gt_row gt_right">259.530</td>
<td headers="Frontal  Anova (p)" class="gt_row gt_right">0.007</td>
<td headers="Frontal  q Value" class="gt_row gt_right">0.024</td>
<td headers="Frontal  Max fold change" class="gt_row gt_right">1.733</td>
<td headers="Frontal  Power" class="gt_row gt_right">0.863</td>
<td headers="Frontal  Highest mean condition" class="gt_row gt_left">2A</td>
<td headers="Frontal  Lowest mean condition" class="gt_row gt_left">1A</td>
<td headers="Frontal  highest_mean" class="gt_row gt_left">POS</td>
<td headers="Frontal  lowest_mean" class="gt_row gt_left">POS</td>
<td headers="Frontal  Mass" class="gt_row gt_right">98252</td>
<td headers="Frontal  Description" class="gt_row gt_left">Dynamin-3 OS=Rattus norvegicus GN=Dnm3 PE=1 SV=2</td></tr>
    <tr><td headers="Frontal  Accession" class="gt_row gt_left">THIL_RAT</td>
<td headers="Frontal  Peptide count" class="gt_row gt_right">9</td>
<td headers="Frontal  Unique peptides" class="gt_row gt_right">1</td>
<td headers="Frontal  Confidence score" class="gt_row gt_right">156.950</td>
<td headers="Frontal  Anova (p)" class="gt_row gt_right">0.012</td>
<td headers="Frontal  q Value" class="gt_row gt_right">0.024</td>
<td headers="Frontal  Max fold change" class="gt_row gt_right">1.922</td>
<td headers="Frontal  Power" class="gt_row gt_right">0.809</td>
<td headers="Frontal  Highest mean condition" class="gt_row gt_left">2A</td>
<td headers="Frontal  Lowest mean condition" class="gt_row gt_left">1A</td>
<td headers="Frontal  highest_mean" class="gt_row gt_left">POS</td>
<td headers="Frontal  lowest_mean" class="gt_row gt_left">POS</td>
<td headers="Frontal  Mass" class="gt_row gt_right">45009</td>
<td headers="Frontal  Description" class="gt_row gt_left">Acetyl-CoA acetyltransferase, mitochondrial OS=Rattus norvegicus GN=Acat1 PE=1 SV=1</td></tr>
    <tr><td headers="Frontal  Accession" class="gt_row gt_left">PGRC1_MOUSE</td>
<td headers="Frontal  Peptide count" class="gt_row gt_right">3</td>
<td headers="Frontal  Unique peptides" class="gt_row gt_right">2</td>
<td headers="Frontal  Confidence score" class="gt_row gt_right">83.180</td>
<td headers="Frontal  Anova (p)" class="gt_row gt_right">0.018</td>
<td headers="Frontal  q Value" class="gt_row gt_right">0.024</td>
<td headers="Frontal  Max fold change" class="gt_row gt_right">2.764</td>
<td headers="Frontal  Power" class="gt_row gt_right">0.745</td>
<td headers="Frontal  Highest mean condition" class="gt_row gt_left">2A</td>
<td headers="Frontal  Lowest mean condition" class="gt_row gt_left">1A</td>
<td headers="Frontal  highest_mean" class="gt_row gt_left">POS</td>
<td headers="Frontal  lowest_mean" class="gt_row gt_left">POS</td>
<td headers="Frontal  Mass" class="gt_row gt_right">21795</td>
<td headers="Frontal  Description" class="gt_row gt_left">Membrane-associated progesterone receptor component 1 OS=Mus musculus GN=Pgrmc1 PE=1 SV=4</td></tr>
    <tr><td headers="Frontal  Accession" class="gt_row gt_left">VATH_MOUSE</td>
<td headers="Frontal  Peptide count" class="gt_row gt_right">9</td>
<td headers="Frontal  Unique peptides" class="gt_row gt_right">2</td>
<td headers="Frontal  Confidence score" class="gt_row gt_right">125.540</td>
<td headers="Frontal  Anova (p)" class="gt_row gt_right">0.020</td>
<td headers="Frontal  q Value" class="gt_row gt_right">0.024</td>
<td headers="Frontal  Max fold change" class="gt_row gt_right">1.743</td>
<td headers="Frontal  Power" class="gt_row gt_right">0.728</td>
<td headers="Frontal  Highest mean condition" class="gt_row gt_left">2A</td>
<td headers="Frontal  Lowest mean condition" class="gt_row gt_left">1A</td>
<td headers="Frontal  highest_mean" class="gt_row gt_left">POS</td>
<td headers="Frontal  lowest_mean" class="gt_row gt_left">POS</td>
<td headers="Frontal  Mass" class="gt_row gt_right">56275</td>
<td headers="Frontal  Description" class="gt_row gt_left">V-type proton ATPase subunit H OS=Mus musculus GN=Atp6v1h PE=1 SV=1</td></tr>
    <tr><td headers="Frontal  Accession" class="gt_row gt_left">ADDB_MOUSE</td>
<td headers="Frontal  Peptide count" class="gt_row gt_right">9</td>
<td headers="Frontal  Unique peptides" class="gt_row gt_right">0</td>
<td headers="Frontal  Confidence score" class="gt_row gt_right">31.660</td>
<td headers="Frontal  Anova (p)" class="gt_row gt_right">0.020</td>
<td headers="Frontal  q Value" class="gt_row gt_right">0.024</td>
<td headers="Frontal  Max fold change" class="gt_row gt_right">1.770</td>
<td headers="Frontal  Power" class="gt_row gt_right">0.723</td>
<td headers="Frontal  Highest mean condition" class="gt_row gt_left">2A</td>
<td headers="Frontal  Lowest mean condition" class="gt_row gt_left">1A</td>
<td headers="Frontal  highest_mean" class="gt_row gt_left">POS</td>
<td headers="Frontal  lowest_mean" class="gt_row gt_left">POS</td>
<td headers="Frontal  Mass" class="gt_row gt_right">81047</td>
<td headers="Frontal  Description" class="gt_row gt_left">Beta-adducin OS=Mus musculus GN=Add2 PE=1 SV=4</td></tr>
    <tr><td headers="Frontal  Accession" class="gt_row gt_left">NDUS6_MOUSE</td>
<td headers="Frontal  Peptide count" class="gt_row gt_right">3</td>
<td headers="Frontal  Unique peptides" class="gt_row gt_right">0</td>
<td headers="Frontal  Confidence score" class="gt_row gt_right">66.230</td>
<td headers="Frontal  Anova (p)" class="gt_row gt_right">0.020</td>
<td headers="Frontal  q Value" class="gt_row gt_right">0.024</td>
<td headers="Frontal  Max fold change" class="gt_row gt_right">2.345</td>
<td headers="Frontal  Power" class="gt_row gt_right">0.733</td>
<td headers="Frontal  Highest mean condition" class="gt_row gt_left">2A</td>
<td headers="Frontal  Lowest mean condition" class="gt_row gt_left">1A</td>
<td headers="Frontal  highest_mean" class="gt_row gt_left">POS</td>
<td headers="Frontal  lowest_mean" class="gt_row gt_left">POS</td>
<td headers="Frontal  Mass" class="gt_row gt_right">13183</td>
<td headers="Frontal  Description" class="gt_row gt_left">NADH dehydrogenase [ubiquinone] iron-sulfur protein 6, mitochondrial OS=Mus musculus GN=Ndufs6 PE=1 SV=2</td></tr>
    <tr><td headers="Frontal  Accession" class="gt_row gt_left">SDHA_MOUSE</td>
<td headers="Frontal  Peptide count" class="gt_row gt_right">7</td>
<td headers="Frontal  Unique peptides" class="gt_row gt_right">0</td>
<td headers="Frontal  Confidence score" class="gt_row gt_right">75.660</td>
<td headers="Frontal  Anova (p)" class="gt_row gt_right">0.022</td>
<td headers="Frontal  q Value" class="gt_row gt_right">0.024</td>
<td headers="Frontal  Max fold change" class="gt_row gt_right">2.186</td>
<td headers="Frontal  Power" class="gt_row gt_right">0.705</td>
<td headers="Frontal  Highest mean condition" class="gt_row gt_left">2A</td>
<td headers="Frontal  Lowest mean condition" class="gt_row gt_left">1A</td>
<td headers="Frontal  highest_mean" class="gt_row gt_left">POS</td>
<td headers="Frontal  lowest_mean" class="gt_row gt_left">POS</td>
<td headers="Frontal  Mass" class="gt_row gt_right">73623</td>
<td headers="Frontal  Description" class="gt_row gt_left">Succinate dehydrogenase [ubiquinone] flavoprotein subunit, mitochondrial OS=Mus musculus GN=Sdha PE=1 SV=1</td></tr>
    <tr><td headers="Frontal  Accession" class="gt_row gt_left">AQP4_RAT;AQP4_MOUSE</td>
<td headers="Frontal  Peptide count" class="gt_row gt_right">3</td>
<td headers="Frontal  Unique peptides" class="gt_row gt_right">3</td>
<td headers="Frontal  Confidence score" class="gt_row gt_right">49.090</td>
<td headers="Frontal  Anova (p)" class="gt_row gt_right">0.023</td>
<td headers="Frontal  q Value" class="gt_row gt_right">0.024</td>
<td headers="Frontal  Max fold change" class="gt_row gt_right">2.098</td>
<td headers="Frontal  Power" class="gt_row gt_right">0.700</td>
<td headers="Frontal  Highest mean condition" class="gt_row gt_left">2A</td>
<td headers="Frontal  Lowest mean condition" class="gt_row gt_left">1A</td>
<td headers="Frontal  highest_mean" class="gt_row gt_left">POS</td>
<td headers="Frontal  lowest_mean" class="gt_row gt_left">POS</td>
<td headers="Frontal  Mass" class="gt_row gt_right">34914</td>
<td headers="Frontal  Description" class="gt_row gt_left">Aquaporin-4 OS=Rattus norvegicus GN=Aqp4 PE=1 SV=1</td></tr>
    <tr><td headers="Frontal  Accession" class="gt_row gt_left">GRIA2_RAT</td>
<td headers="Frontal  Peptide count" class="gt_row gt_right">12</td>
<td headers="Frontal  Unique peptides" class="gt_row gt_right">3</td>
<td headers="Frontal  Confidence score" class="gt_row gt_right">69.790</td>
<td headers="Frontal  Anova (p)" class="gt_row gt_right">0.030</td>
<td headers="Frontal  q Value" class="gt_row gt_right">0.024</td>
<td headers="Frontal  Max fold change" class="gt_row gt_right">2.617</td>
<td headers="Frontal  Power" class="gt_row gt_right">0.672</td>
<td headers="Frontal  Highest mean condition" class="gt_row gt_left">2A</td>
<td headers="Frontal  Lowest mean condition" class="gt_row gt_left">1A</td>
<td headers="Frontal  highest_mean" class="gt_row gt_left">POS</td>
<td headers="Frontal  lowest_mean" class="gt_row gt_left">POS</td>
<td headers="Frontal  Mass" class="gt_row gt_right">99252</td>
<td headers="Frontal  Description" class="gt_row gt_left">Glutamate receptor 2 OS=Rattus norvegicus GN=Gria2 PE=1 SV=2</td></tr>
    <tr><td headers="Frontal  Accession" class="gt_row gt_left">STX2_RAT</td>
<td headers="Frontal  Peptide count" class="gt_row gt_right">2</td>
<td headers="Frontal  Unique peptides" class="gt_row gt_right">1</td>
<td headers="Frontal  Confidence score" class="gt_row gt_right">43.120</td>
<td headers="Frontal  Anova (p)" class="gt_row gt_right">0.030</td>
<td headers="Frontal  q Value" class="gt_row gt_right">0.024</td>
<td headers="Frontal  Max fold change" class="gt_row gt_right">4.534</td>
<td headers="Frontal  Power" class="gt_row gt_right">0.650</td>
<td headers="Frontal  Highest mean condition" class="gt_row gt_left">2A</td>
<td headers="Frontal  Lowest mean condition" class="gt_row gt_left">1A</td>
<td headers="Frontal  highest_mean" class="gt_row gt_left">POS</td>
<td headers="Frontal  lowest_mean" class="gt_row gt_left">POS</td>
<td headers="Frontal  Mass" class="gt_row gt_right">33452</td>
<td headers="Frontal  Description" class="gt_row gt_left">Syntaxin-2 OS=Rattus norvegicus GN=Stx2 PE=1 SV=2</td></tr>
    <tr><td headers="Frontal  Accession" class="gt_row gt_left">OMP_RAT</td>
<td headers="Frontal  Peptide count" class="gt_row gt_right">4</td>
<td headers="Frontal  Unique peptides" class="gt_row gt_right">2</td>
<td headers="Frontal  Confidence score" class="gt_row gt_right">93.210</td>
<td headers="Frontal  Anova (p)" class="gt_row gt_right">0.032</td>
<td headers="Frontal  q Value" class="gt_row gt_right">0.024</td>
<td headers="Frontal  Max fold change" class="gt_row gt_right">3.617</td>
<td headers="Frontal  Power" class="gt_row gt_right">0.641</td>
<td headers="Frontal  Highest mean condition" class="gt_row gt_left">2A</td>
<td headers="Frontal  Lowest mean condition" class="gt_row gt_left">1A</td>
<td headers="Frontal  highest_mean" class="gt_row gt_left">POS</td>
<td headers="Frontal  lowest_mean" class="gt_row gt_left">POS</td>
<td headers="Frontal  Mass" class="gt_row gt_right">18841</td>
<td headers="Frontal  Description" class="gt_row gt_left">Olfactory marker protein OS=Rattus norvegicus GN=Omp PE=1 SV=2</td></tr>
    <tr><td headers="Frontal  Accession" class="gt_row gt_left">ROA1_MOUSE</td>
<td headers="Frontal  Peptide count" class="gt_row gt_right">6</td>
<td headers="Frontal  Unique peptides" class="gt_row gt_right">3</td>
<td headers="Frontal  Confidence score" class="gt_row gt_right">158.280</td>
<td headers="Frontal  Anova (p)" class="gt_row gt_right">0.032</td>
<td headers="Frontal  q Value" class="gt_row gt_right">0.024</td>
<td headers="Frontal  Max fold change" class="gt_row gt_right">2.036</td>
<td headers="Frontal  Power" class="gt_row gt_right">0.643</td>
<td headers="Frontal  Highest mean condition" class="gt_row gt_left">2A</td>
<td headers="Frontal  Lowest mean condition" class="gt_row gt_left">1A</td>
<td headers="Frontal  highest_mean" class="gt_row gt_left">POS</td>
<td headers="Frontal  lowest_mean" class="gt_row gt_left">POS</td>
<td headers="Frontal  Mass" class="gt_row gt_right">34289</td>
<td headers="Frontal  Description" class="gt_row gt_left">Heterogeneous nuclear ribonucleoprotein A1 OS=Mus musculus GN=Hnrnpa1 PE=1 SV=2</td></tr>
    <tr><td headers="Frontal  Accession" class="gt_row gt_left">HNRPD_RAT</td>
<td headers="Frontal  Peptide count" class="gt_row gt_right">8</td>
<td headers="Frontal  Unique peptides" class="gt_row gt_right">3</td>
<td headers="Frontal  Confidence score" class="gt_row gt_right">100.830</td>
<td headers="Frontal  Anova (p)" class="gt_row gt_right">0.033</td>
<td headers="Frontal  q Value" class="gt_row gt_right">0.024</td>
<td headers="Frontal  Max fold change" class="gt_row gt_right">2.410</td>
<td headers="Frontal  Power" class="gt_row gt_right">0.647</td>
<td headers="Frontal  Highest mean condition" class="gt_row gt_left">2A</td>
<td headers="Frontal  Lowest mean condition" class="gt_row gt_left">1A</td>
<td headers="Frontal  highest_mean" class="gt_row gt_left">POS</td>
<td headers="Frontal  lowest_mean" class="gt_row gt_left">POS</td>
<td headers="Frontal  Mass" class="gt_row gt_right">38365</td>
<td headers="Frontal  Description" class="gt_row gt_left">Heterogeneous nuclear ribonucleoprotein D0 OS=Rattus norvegicus GN=Hnrnpd PE=1 SV=2</td></tr>
    <tr><td headers="Frontal  Accession" class="gt_row gt_left">BD1L1_MOUSE</td>
<td headers="Frontal  Peptide count" class="gt_row gt_right">32</td>
<td headers="Frontal  Unique peptides" class="gt_row gt_right">5</td>
<td headers="Frontal  Confidence score" class="gt_row gt_right">251.290</td>
<td headers="Frontal  Anova (p)" class="gt_row gt_right">0.033</td>
<td headers="Frontal  q Value" class="gt_row gt_right">0.024</td>
<td headers="Frontal  Max fold change" class="gt_row gt_right">2.664</td>
<td headers="Frontal  Power" class="gt_row gt_right">0.649</td>
<td headers="Frontal  Highest mean condition" class="gt_row gt_left">2A</td>
<td headers="Frontal  Lowest mean condition" class="gt_row gt_left">1A</td>
<td headers="Frontal  highest_mean" class="gt_row gt_left">POS</td>
<td headers="Frontal  lowest_mean" class="gt_row gt_left">POS</td>
<td headers="Frontal  Mass" class="gt_row gt_right">330220</td>
<td headers="Frontal  Description" class="gt_row gt_left">Biorientation of chromosomes in cell division protein 1-like 1 OS=Mus musculus GN=Bod1l PE=1 SV=1</td></tr>
    <tr><td headers="Frontal  Accession" class="gt_row gt_left">VATD_MOUSE</td>
<td headers="Frontal  Peptide count" class="gt_row gt_right">7</td>
<td headers="Frontal  Unique peptides" class="gt_row gt_right">3</td>
<td headers="Frontal  Confidence score" class="gt_row gt_right">65.340</td>
<td headers="Frontal  Anova (p)" class="gt_row gt_right">0.033</td>
<td headers="Frontal  q Value" class="gt_row gt_right">0.024</td>
<td headers="Frontal  Max fold change" class="gt_row gt_right">2.365</td>
<td headers="Frontal  Power" class="gt_row gt_right">0.689</td>
<td headers="Frontal  Highest mean condition" class="gt_row gt_left">2A</td>
<td headers="Frontal  Lowest mean condition" class="gt_row gt_left">3A</td>
<td headers="Frontal  highest_mean" class="gt_row gt_left">POS</td>
<td headers="Frontal  lowest_mean" class="gt_row gt_left">POS</td>
<td headers="Frontal  Mass" class="gt_row gt_right">28351</td>
<td headers="Frontal  Description" class="gt_row gt_left">V-type proton ATPase subunit D OS=Mus musculus GN=Atp6v1d PE=1 SV=1</td></tr>
    <tr><td headers="Frontal  Accession" class="gt_row gt_left">STAG3_MOUSE</td>
<td headers="Frontal  Peptide count" class="gt_row gt_right">7</td>
<td headers="Frontal  Unique peptides" class="gt_row gt_right">1</td>
<td headers="Frontal  Confidence score" class="gt_row gt_right">30.690</td>
<td headers="Frontal  Anova (p)" class="gt_row gt_right">0.033</td>
<td headers="Frontal  q Value" class="gt_row gt_right">0.024</td>
<td headers="Frontal  Max fold change" class="gt_row gt_right">1.809</td>
<td headers="Frontal  Power" class="gt_row gt_right">0.635</td>
<td headers="Frontal  Highest mean condition" class="gt_row gt_left">3A</td>
<td headers="Frontal  Lowest mean condition" class="gt_row gt_left">1A</td>
<td headers="Frontal  highest_mean" class="gt_row gt_left">POS</td>
<td headers="Frontal  lowest_mean" class="gt_row gt_left">POS</td>
<td headers="Frontal  Mass" class="gt_row gt_right">142731</td>
<td headers="Frontal  Description" class="gt_row gt_left">Cohesin subunit SA-3 OS=Mus musculus GN=Stag3 PE=1 SV=2</td></tr>
    <tr><td headers="Frontal  Accession" class="gt_row gt_left">CAPZB_RAT</td>
<td headers="Frontal  Peptide count" class="gt_row gt_right">4</td>
<td headers="Frontal  Unique peptides" class="gt_row gt_right">1</td>
<td headers="Frontal  Confidence score" class="gt_row gt_right">71.340</td>
<td headers="Frontal  Anova (p)" class="gt_row gt_right">0.035</td>
<td headers="Frontal  q Value" class="gt_row gt_right">0.024</td>
<td headers="Frontal  Max fold change" class="gt_row gt_right">1.584</td>
<td headers="Frontal  Power" class="gt_row gt_right">0.625</td>
<td headers="Frontal  Highest mean condition" class="gt_row gt_left">2A</td>
<td headers="Frontal  Lowest mean condition" class="gt_row gt_left">1A</td>
<td headers="Frontal  highest_mean" class="gt_row gt_left">POS</td>
<td headers="Frontal  lowest_mean" class="gt_row gt_left">POS</td>
<td headers="Frontal  Mass" class="gt_row gt_right">30952</td>
<td headers="Frontal  Description" class="gt_row gt_left">F-actin-capping protein subunit beta OS=Rattus norvegicus GN=Capzb PE=1 SV=1</td></tr>
    <tr><td headers="Frontal  Accession" class="gt_row gt_left">BACH_MOUSE</td>
<td headers="Frontal  Peptide count" class="gt_row gt_right">5</td>
<td headers="Frontal  Unique peptides" class="gt_row gt_right">1</td>
<td headers="Frontal  Confidence score" class="gt_row gt_right">96.150</td>
<td headers="Frontal  Anova (p)" class="gt_row gt_right">0.036</td>
<td headers="Frontal  q Value" class="gt_row gt_right">0.024</td>
<td headers="Frontal  Max fold change" class="gt_row gt_right">4.021</td>
<td headers="Frontal  Power" class="gt_row gt_right">0.637</td>
<td headers="Frontal  Highest mean condition" class="gt_row gt_left">2A</td>
<td headers="Frontal  Lowest mean condition" class="gt_row gt_left">1A</td>
<td headers="Frontal  highest_mean" class="gt_row gt_left">POS</td>
<td headers="Frontal  lowest_mean" class="gt_row gt_left">POS</td>
<td headers="Frontal  Mass" class="gt_row gt_right">42966</td>
<td headers="Frontal  Description" class="gt_row gt_left">Cytosolic acyl coenzyme A thioester hydrolase OS=Mus musculus GN=Acot7 PE=1 SV=2</td></tr>
    <tr><td headers="Frontal  Accession" class="gt_row gt_left">MOG_MOUSE</td>
<td headers="Frontal  Peptide count" class="gt_row gt_right">1</td>
<td headers="Frontal  Unique peptides" class="gt_row gt_right">1</td>
<td headers="Frontal  Confidence score" class="gt_row gt_right">66.330</td>
<td headers="Frontal  Anova (p)" class="gt_row gt_right">0.037</td>
<td headers="Frontal  q Value" class="gt_row gt_right">0.024</td>
<td headers="Frontal  Max fold change" class="gt_row gt_right">3.050</td>
<td headers="Frontal  Power" class="gt_row gt_right">0.631</td>
<td headers="Frontal  Highest mean condition" class="gt_row gt_left">2A</td>
<td headers="Frontal  Lowest mean condition" class="gt_row gt_left">1A</td>
<td headers="Frontal  highest_mean" class="gt_row gt_left">POS</td>
<td headers="Frontal  lowest_mean" class="gt_row gt_left">POS</td>
<td headers="Frontal  Mass" class="gt_row gt_right">28652</td>
<td headers="Frontal  Description" class="gt_row gt_left">Myelin-oligodendrocyte glycoprotein OS=Mus musculus GN=Mog PE=1 SV=1</td></tr>
    <tr><td headers="Frontal  Accession" class="gt_row gt_left">HCDH_RAT</td>
<td headers="Frontal  Peptide count" class="gt_row gt_right">4</td>
<td headers="Frontal  Unique peptides" class="gt_row gt_right">0</td>
<td headers="Frontal  Confidence score" class="gt_row gt_right">68.350</td>
<td headers="Frontal  Anova (p)" class="gt_row gt_right">0.038</td>
<td headers="Frontal  q Value" class="gt_row gt_right">0.024</td>
<td headers="Frontal  Max fold change" class="gt_row gt_right">2.626</td>
<td headers="Frontal  Power" class="gt_row gt_right">0.612</td>
<td headers="Frontal  Highest mean condition" class="gt_row gt_left">2A</td>
<td headers="Frontal  Lowest mean condition" class="gt_row gt_left">1A</td>
<td headers="Frontal  highest_mean" class="gt_row gt_left">POS</td>
<td headers="Frontal  lowest_mean" class="gt_row gt_left">POS</td>
<td headers="Frontal  Mass" class="gt_row gt_right">34540</td>
<td headers="Frontal  Description" class="gt_row gt_left">Hydroxyacyl-coenzyme A dehydrogenase, mitochondrial OS=Rattus norvegicus GN=Hadh PE=2 SV=1</td></tr>
    <tr><td headers="Frontal  Accession" class="gt_row gt_left">SEP11_MOUSE;SEP10_RAT;SEPT8_RAT</td>
<td headers="Frontal  Peptide count" class="gt_row gt_right">4</td>
<td headers="Frontal  Unique peptides" class="gt_row gt_right">0</td>
<td headers="Frontal  Confidence score" class="gt_row gt_right">81.400</td>
<td headers="Frontal  Anova (p)" class="gt_row gt_right">0.039</td>
<td headers="Frontal  q Value" class="gt_row gt_right">0.024</td>
<td headers="Frontal  Max fold change" class="gt_row gt_right">1.937</td>
<td headers="Frontal  Power" class="gt_row gt_right">0.624</td>
<td headers="Frontal  Highest mean condition" class="gt_row gt_left">2A</td>
<td headers="Frontal  Lowest mean condition" class="gt_row gt_left">1A</td>
<td headers="Frontal  highest_mean" class="gt_row gt_left">POS</td>
<td headers="Frontal  lowest_mean" class="gt_row gt_left">POS</td>
<td headers="Frontal  Mass" class="gt_row gt_right">50005</td>
<td headers="Frontal  Description" class="gt_row gt_left">Septin-11 OS=Mus musculus GN=Sept11 PE=1 SV=4</td></tr>
    <tr><td headers="Frontal  Accession" class="gt_row gt_left">KPCE_MOUSE</td>
<td headers="Frontal  Peptide count" class="gt_row gt_right">5</td>
<td headers="Frontal  Unique peptides" class="gt_row gt_right">1</td>
<td headers="Frontal  Confidence score" class="gt_row gt_right">79.300</td>
<td headers="Frontal  Anova (p)" class="gt_row gt_right">0.040</td>
<td headers="Frontal  q Value" class="gt_row gt_right">0.024</td>
<td headers="Frontal  Max fold change" class="gt_row gt_right">2.009</td>
<td headers="Frontal  Power" class="gt_row gt_right">0.630</td>
<td headers="Frontal  Highest mean condition" class="gt_row gt_left">2A</td>
<td headers="Frontal  Lowest mean condition" class="gt_row gt_left">1A</td>
<td headers="Frontal  highest_mean" class="gt_row gt_left">POS</td>
<td headers="Frontal  lowest_mean" class="gt_row gt_left">POS</td>
<td headers="Frontal  Mass" class="gt_row gt_right">84875</td>
<td headers="Frontal  Description" class="gt_row gt_left">Protein kinase C epsilon type OS=Mus musculus GN=Prkce PE=1 SV=1</td></tr>
    <tr><td headers="Frontal  Accession" class="gt_row gt_left">QCR2_RAT</td>
<td headers="Frontal  Peptide count" class="gt_row gt_right">15</td>
<td headers="Frontal  Unique peptides" class="gt_row gt_right">0</td>
<td headers="Frontal  Confidence score" class="gt_row gt_right">477.820</td>
<td headers="Frontal  Anova (p)" class="gt_row gt_right">0.040</td>
<td headers="Frontal  q Value" class="gt_row gt_right">0.024</td>
<td headers="Frontal  Max fold change" class="gt_row gt_right">2.750</td>
<td headers="Frontal  Power" class="gt_row gt_right">0.622</td>
<td headers="Frontal  Highest mean condition" class="gt_row gt_left">2A</td>
<td headers="Frontal  Lowest mean condition" class="gt_row gt_left">1A</td>
<td headers="Frontal  highest_mean" class="gt_row gt_left">POS</td>
<td headers="Frontal  lowest_mean" class="gt_row gt_left">POS</td>
<td headers="Frontal  Mass" class="gt_row gt_right">48423</td>
<td headers="Frontal  Description" class="gt_row gt_left">Cytochrome b-c1 complex subunit 2, mitochondrial OS=Rattus norvegicus GN=Uqcrc2 PE=1 SV=2</td></tr>
    <tr><td headers="Frontal  Accession" class="gt_row gt_left">INP4A_RAT;INP4A_MOUSE</td>
<td headers="Frontal  Peptide count" class="gt_row gt_right">8</td>
<td headers="Frontal  Unique peptides" class="gt_row gt_right">2</td>
<td headers="Frontal  Confidence score" class="gt_row gt_right">68.830</td>
<td headers="Frontal  Anova (p)" class="gt_row gt_right">0.041</td>
<td headers="Frontal  q Value" class="gt_row gt_right">0.024</td>
<td headers="Frontal  Max fold change" class="gt_row gt_right">1.977</td>
<td headers="Frontal  Power" class="gt_row gt_right">0.624</td>
<td headers="Frontal  Highest mean condition" class="gt_row gt_left">2A</td>
<td headers="Frontal  Lowest mean condition" class="gt_row gt_left">1A</td>
<td headers="Frontal  highest_mean" class="gt_row gt_left">POS</td>
<td headers="Frontal  lowest_mean" class="gt_row gt_left">POS</td>
<td headers="Frontal  Mass" class="gt_row gt_right">107004</td>
<td headers="Frontal  Description" class="gt_row gt_left">Type I inositol 3,4-bisphosphate 4-phosphatase OS=Rattus norvegicus GN=Inpp4a PE=1 SV=1</td></tr>
    <tr><td headers="Frontal  Accession" class="gt_row gt_left">CMC1_MOUSE</td>
<td headers="Frontal  Peptide count" class="gt_row gt_right">10</td>
<td headers="Frontal  Unique peptides" class="gt_row gt_right">3</td>
<td headers="Frontal  Confidence score" class="gt_row gt_right">383.050</td>
<td headers="Frontal  Anova (p)" class="gt_row gt_right">0.041</td>
<td headers="Frontal  q Value" class="gt_row gt_right">0.024</td>
<td headers="Frontal  Max fold change" class="gt_row gt_right">1.793</td>
<td headers="Frontal  Power" class="gt_row gt_right">0.607</td>
<td headers="Frontal  Highest mean condition" class="gt_row gt_left">2A</td>
<td headers="Frontal  Lowest mean condition" class="gt_row gt_left">1A</td>
<td headers="Frontal  highest_mean" class="gt_row gt_left">POS</td>
<td headers="Frontal  lowest_mean" class="gt_row gt_left">POS</td>
<td headers="Frontal  Mass" class="gt_row gt_right">74922</td>
<td headers="Frontal  Description" class="gt_row gt_left">Calcium-binding mitochondrial carrier protein Aralar1 OS=Mus musculus GN=Slc25a12 PE=1 SV=1</td></tr>
    <tr><td headers="Frontal  Accession" class="gt_row gt_left">MYPR_MOUSE</td>
<td headers="Frontal  Peptide count" class="gt_row gt_right">10</td>
<td headers="Frontal  Unique peptides" class="gt_row gt_right">1</td>
<td headers="Frontal  Confidence score" class="gt_row gt_right">370.010</td>
<td headers="Frontal  Anova (p)" class="gt_row gt_right">0.041</td>
<td headers="Frontal  q Value" class="gt_row gt_right">0.024</td>
<td headers="Frontal  Max fold change" class="gt_row gt_right">2.617</td>
<td headers="Frontal  Power" class="gt_row gt_right">0.597</td>
<td headers="Frontal  Highest mean condition" class="gt_row gt_left">2A</td>
<td headers="Frontal  Lowest mean condition" class="gt_row gt_left">1A</td>
<td headers="Frontal  highest_mean" class="gt_row gt_left">POS</td>
<td headers="Frontal  lowest_mean" class="gt_row gt_left">POS</td>
<td headers="Frontal  Mass" class="gt_row gt_right">30855</td>
<td headers="Frontal  Description" class="gt_row gt_left">Myelin proteolipid protein OS=Mus musculus GN=Plp1 PE=1 SV=2</td></tr>
    <tr><td headers="Frontal  Accession" class="gt_row gt_left">SPA3M_MOUSE</td>
<td headers="Frontal  Peptide count" class="gt_row gt_right">3</td>
<td headers="Frontal  Unique peptides" class="gt_row gt_right">1</td>
<td headers="Frontal  Confidence score" class="gt_row gt_right">214.420</td>
<td headers="Frontal  Anova (p)" class="gt_row gt_right">0.041</td>
<td headers="Frontal  q Value" class="gt_row gt_right">0.024</td>
<td headers="Frontal  Max fold change" class="gt_row gt_right">2.583</td>
<td headers="Frontal  Power" class="gt_row gt_right">0.599</td>
<td headers="Frontal  Highest mean condition" class="gt_row gt_left">2A</td>
<td headers="Frontal  Lowest mean condition" class="gt_row gt_left">1A</td>
<td headers="Frontal  highest_mean" class="gt_row gt_left">POS</td>
<td headers="Frontal  lowest_mean" class="gt_row gt_left">POS</td>
<td headers="Frontal  Mass" class="gt_row gt_right">47205</td>
<td headers="Frontal  Description" class="gt_row gt_left">Serine protease inhibitor A3M OS=Mus musculus GN=Serpina3m PE=1 SV=2</td></tr>
    <tr><td headers="Frontal  Accession" class="gt_row gt_left">ACTN3_MOUSE</td>
<td headers="Frontal  Peptide count" class="gt_row gt_right">11</td>
<td headers="Frontal  Unique peptides" class="gt_row gt_right">2</td>
<td headers="Frontal  Confidence score" class="gt_row gt_right">120.310</td>
<td headers="Frontal  Anova (p)" class="gt_row gt_right">0.041</td>
<td headers="Frontal  q Value" class="gt_row gt_right">0.024</td>
<td headers="Frontal  Max fold change" class="gt_row gt_right">4.097</td>
<td headers="Frontal  Power" class="gt_row gt_right">0.612</td>
<td headers="Frontal  Highest mean condition" class="gt_row gt_left">2A</td>
<td headers="Frontal  Lowest mean condition" class="gt_row gt_left">1A</td>
<td headers="Frontal  highest_mean" class="gt_row gt_left">POS</td>
<td headers="Frontal  lowest_mean" class="gt_row gt_left">POS</td>
<td headers="Frontal  Mass" class="gt_row gt_right">103605</td>
<td headers="Frontal  Description" class="gt_row gt_left">Alpha-actinin-3 OS=Mus musculus GN=Actn3 PE=2 SV=1</td></tr>
    <tr><td headers="Frontal  Accession" class="gt_row gt_left">NDUV1_MOUSE</td>
<td headers="Frontal  Peptide count" class="gt_row gt_right">6</td>
<td headers="Frontal  Unique peptides" class="gt_row gt_right">2</td>
<td headers="Frontal  Confidence score" class="gt_row gt_right">76.190</td>
<td headers="Frontal  Anova (p)" class="gt_row gt_right">0.042</td>
<td headers="Frontal  q Value" class="gt_row gt_right">0.024</td>
<td headers="Frontal  Max fold change" class="gt_row gt_right">4.163</td>
<td headers="Frontal  Power" class="gt_row gt_right">0.604</td>
<td headers="Frontal  Highest mean condition" class="gt_row gt_left">2A</td>
<td headers="Frontal  Lowest mean condition" class="gt_row gt_left">1A</td>
<td headers="Frontal  highest_mean" class="gt_row gt_left">POS</td>
<td headers="Frontal  lowest_mean" class="gt_row gt_left">POS</td>
<td headers="Frontal  Mass" class="gt_row gt_right">51486</td>
<td headers="Frontal  Description" class="gt_row gt_left">NADH dehydrogenase [ubiquinone] flavoprotein 1, mitochondrial OS=Mus musculus GN=Ndufv1 PE=1 SV=1</td></tr>
    <tr><td headers="Frontal  Accession" class="gt_row gt_left">ILF3_MOUSE</td>
<td headers="Frontal  Peptide count" class="gt_row gt_right">9</td>
<td headers="Frontal  Unique peptides" class="gt_row gt_right">1</td>
<td headers="Frontal  Confidence score" class="gt_row gt_right">37.410</td>
<td headers="Frontal  Anova (p)" class="gt_row gt_right">0.043</td>
<td headers="Frontal  q Value" class="gt_row gt_right">0.024</td>
<td headers="Frontal  Max fold change" class="gt_row gt_right">2.288</td>
<td headers="Frontal  Power" class="gt_row gt_right">0.597</td>
<td headers="Frontal  Highest mean condition" class="gt_row gt_left">2A</td>
<td headers="Frontal  Lowest mean condition" class="gt_row gt_left">1A</td>
<td headers="Frontal  highest_mean" class="gt_row gt_left">POS</td>
<td headers="Frontal  lowest_mean" class="gt_row gt_left">POS</td>
<td headers="Frontal  Mass" class="gt_row gt_right">96360</td>
<td headers="Frontal  Description" class="gt_row gt_left">Interleukin enhancer-binding factor 3 OS=Mus musculus GN=Ilf3 PE=1 SV=2</td></tr>
    <tr><td headers="Frontal  Accession" class="gt_row gt_left">VATF_MOUSE</td>
<td headers="Frontal  Peptide count" class="gt_row gt_right">4</td>
<td headers="Frontal  Unique peptides" class="gt_row gt_right">1</td>
<td headers="Frontal  Confidence score" class="gt_row gt_right">220.520</td>
<td headers="Frontal  Anova (p)" class="gt_row gt_right">0.044</td>
<td headers="Frontal  q Value" class="gt_row gt_right">0.024</td>
<td headers="Frontal  Max fold change" class="gt_row gt_right">1.878</td>
<td headers="Frontal  Power" class="gt_row gt_right">0.620</td>
<td headers="Frontal  Highest mean condition" class="gt_row gt_left">2A</td>
<td headers="Frontal  Lowest mean condition" class="gt_row gt_left">3A</td>
<td headers="Frontal  highest_mean" class="gt_row gt_left">POS</td>
<td headers="Frontal  lowest_mean" class="gt_row gt_left">POS</td>
<td headers="Frontal  Mass" class="gt_row gt_right">13362</td>
<td headers="Frontal  Description" class="gt_row gt_left">V-type proton ATPase subunit F OS=Mus musculus GN=Atp6v1f PE=1 SV=2</td></tr>
    <tr><td headers="Frontal  Accession" class="gt_row gt_left">CISY_MOUSE</td>
<td headers="Frontal  Peptide count" class="gt_row gt_right">13</td>
<td headers="Frontal  Unique peptides" class="gt_row gt_right">3</td>
<td headers="Frontal  Confidence score" class="gt_row gt_right">399.480</td>
<td headers="Frontal  Anova (p)" class="gt_row gt_right">0.045</td>
<td headers="Frontal  q Value" class="gt_row gt_right">0.024</td>
<td headers="Frontal  Max fold change" class="gt_row gt_right">2.011</td>
<td headers="Frontal  Power" class="gt_row gt_right">0.583</td>
<td headers="Frontal  Highest mean condition" class="gt_row gt_left">2A</td>
<td headers="Frontal  Lowest mean condition" class="gt_row gt_left">1A</td>
<td headers="Frontal  highest_mean" class="gt_row gt_left">POS</td>
<td headers="Frontal  lowest_mean" class="gt_row gt_left">POS</td>
<td headers="Frontal  Mass" class="gt_row gt_right">51988</td>
<td headers="Frontal  Description" class="gt_row gt_left">Citrate synthase, mitochondrial OS=Mus musculus GN=Cs PE=1 SV=1</td></tr>
    <tr><td headers="Frontal  Accession" class="gt_row gt_left">BACH_RAT</td>
<td headers="Frontal  Peptide count" class="gt_row gt_right">3</td>
<td headers="Frontal  Unique peptides" class="gt_row gt_right">0</td>
<td headers="Frontal  Confidence score" class="gt_row gt_right">84.570</td>
<td headers="Frontal  Anova (p)" class="gt_row gt_right">0.045</td>
<td headers="Frontal  q Value" class="gt_row gt_right">0.024</td>
<td headers="Frontal  Max fold change" class="gt_row gt_right">1.712</td>
<td headers="Frontal  Power" class="gt_row gt_right">0.580</td>
<td headers="Frontal  Highest mean condition" class="gt_row gt_left">3A</td>
<td headers="Frontal  Lowest mean condition" class="gt_row gt_left">1A</td>
<td headers="Frontal  highest_mean" class="gt_row gt_left">POS</td>
<td headers="Frontal  lowest_mean" class="gt_row gt_left">POS</td>
<td headers="Frontal  Mass" class="gt_row gt_right">43164</td>
<td headers="Frontal  Description" class="gt_row gt_left">Cytosolic acyl coenzyme A thioester hydrolase OS=Rattus norvegicus GN=Acot7 PE=1 SV=4</td></tr>
    <tr><td headers="Frontal  Accession" class="gt_row gt_left">H2AX_MOUSE</td>
<td headers="Frontal  Peptide count" class="gt_row gt_right">7</td>
<td headers="Frontal  Unique peptides" class="gt_row gt_right">1</td>
<td headers="Frontal  Confidence score" class="gt_row gt_right">641.100</td>
<td headers="Frontal  Anova (p)" class="gt_row gt_right">0.046</td>
<td headers="Frontal  q Value" class="gt_row gt_right">0.024</td>
<td headers="Frontal  Max fold change" class="gt_row gt_right">2.459</td>
<td headers="Frontal  Power" class="gt_row gt_right">0.585</td>
<td headers="Frontal  Highest mean condition" class="gt_row gt_left">2A</td>
<td headers="Frontal  Lowest mean condition" class="gt_row gt_left">1A</td>
<td headers="Frontal  highest_mean" class="gt_row gt_left">POS</td>
<td headers="Frontal  lowest_mean" class="gt_row gt_left">POS</td>
<td headers="Frontal  Mass" class="gt_row gt_right">15133</td>
<td headers="Frontal  Description" class="gt_row gt_left">Histone H2AX OS=Mus musculus GN=H2afx PE=1 SV=2</td></tr>
    <tr><td headers="Frontal  Accession" class="gt_row gt_left">C1QBP_MOUSE</td>
<td headers="Frontal  Peptide count" class="gt_row gt_right">1</td>
<td headers="Frontal  Unique peptides" class="gt_row gt_right">0</td>
<td headers="Frontal  Confidence score" class="gt_row gt_right">25.090</td>
<td headers="Frontal  Anova (p)" class="gt_row gt_right">0.046</td>
<td headers="Frontal  q Value" class="gt_row gt_right">0.024</td>
<td headers="Frontal  Max fold change" class="gt_row gt_right">2.739</td>
<td headers="Frontal  Power" class="gt_row gt_right">0.578</td>
<td headers="Frontal  Highest mean condition" class="gt_row gt_left">2A</td>
<td headers="Frontal  Lowest mean condition" class="gt_row gt_left">1A</td>
<td headers="Frontal  highest_mean" class="gt_row gt_left">POS</td>
<td headers="Frontal  lowest_mean" class="gt_row gt_left">POS</td>
<td headers="Frontal  Mass" class="gt_row gt_right">31336</td>
<td headers="Frontal  Description" class="gt_row gt_left">Complement component 1 Q subcomponent-binding protein, mitochondrial OS=Mus musculus GN=C1qbp PE=1 SV=1</td></tr>
    <tr><td headers="Frontal  Accession" class="gt_row gt_left">H31_MOUSE;H3C_MOUSE</td>
<td headers="Frontal  Peptide count" class="gt_row gt_right">1</td>
<td headers="Frontal  Unique peptides" class="gt_row gt_right">0</td>
<td headers="Frontal  Confidence score" class="gt_row gt_right">38.290</td>
<td headers="Frontal  Anova (p)" class="gt_row gt_right">0.046</td>
<td headers="Frontal  q Value" class="gt_row gt_right">0.024</td>
<td headers="Frontal  Max fold change" class="gt_row gt_right">2.739</td>
<td headers="Frontal  Power" class="gt_row gt_right">0.578</td>
<td headers="Frontal  Highest mean condition" class="gt_row gt_left">2A</td>
<td headers="Frontal  Lowest mean condition" class="gt_row gt_left">1A</td>
<td headers="Frontal  highest_mean" class="gt_row gt_left">POS</td>
<td headers="Frontal  lowest_mean" class="gt_row gt_left">POS</td>
<td headers="Frontal  Mass" class="gt_row gt_right">15509</td>
<td headers="Frontal  Description" class="gt_row gt_left">Histone H3.1 OS=Mus musculus GN=Hist1h3a PE=1 SV=2</td></tr>
    <tr><td headers="Frontal  Accession" class="gt_row gt_left">CLCB_MOUSE</td>
<td headers="Frontal  Peptide count" class="gt_row gt_right">3</td>
<td headers="Frontal  Unique peptides" class="gt_row gt_right">2</td>
<td headers="Frontal  Confidence score" class="gt_row gt_right">80.990</td>
<td headers="Frontal  Anova (p)" class="gt_row gt_right">0.047</td>
<td headers="Frontal  q Value" class="gt_row gt_right">0.024</td>
<td headers="Frontal  Max fold change" class="gt_row gt_right">1.995</td>
<td headers="Frontal  Power" class="gt_row gt_right">0.610</td>
<td headers="Frontal  Highest mean condition" class="gt_row gt_left">2A</td>
<td headers="Frontal  Lowest mean condition" class="gt_row gt_left">1A</td>
<td headers="Frontal  highest_mean" class="gt_row gt_left">POS</td>
<td headers="Frontal  lowest_mean" class="gt_row gt_left">POS</td>
<td headers="Frontal  Mass" class="gt_row gt_right">25270</td>
<td headers="Frontal  Description" class="gt_row gt_left">Clathrin light chain B OS=Mus musculus GN=Cltb PE=1 SV=1</td></tr>
    <tr><td headers="Frontal  Accession" class="gt_row gt_left">HCDH_MOUSE</td>
<td headers="Frontal  Peptide count" class="gt_row gt_right">7</td>
<td headers="Frontal  Unique peptides" class="gt_row gt_right">3</td>
<td headers="Frontal  Confidence score" class="gt_row gt_right">79.580</td>
<td headers="Frontal  Anova (p)" class="gt_row gt_right">0.047</td>
<td headers="Frontal  q Value" class="gt_row gt_right">0.024</td>
<td headers="Frontal  Max fold change" class="gt_row gt_right">2.498</td>
<td headers="Frontal  Power" class="gt_row gt_right">0.573</td>
<td headers="Frontal  Highest mean condition" class="gt_row gt_left">2A</td>
<td headers="Frontal  Lowest mean condition" class="gt_row gt_left">1A</td>
<td headers="Frontal  highest_mean" class="gt_row gt_left">POS</td>
<td headers="Frontal  lowest_mean" class="gt_row gt_left">POS</td>
<td headers="Frontal  Mass" class="gt_row gt_right">34613</td>
<td headers="Frontal  Description" class="gt_row gt_left">Hydroxyacyl-coenzyme A dehydrogenase, mitochondrial OS=Mus musculus GN=Hadh PE=1 SV=2</td></tr>
    <tr><td headers="Frontal  Accession" class="gt_row gt_left">NAC2_RAT</td>
<td headers="Frontal  Peptide count" class="gt_row gt_right">10</td>
<td headers="Frontal  Unique peptides" class="gt_row gt_right">1</td>
<td headers="Frontal  Confidence score" class="gt_row gt_right">133.140</td>
<td headers="Frontal  Anova (p)" class="gt_row gt_right">0.049</td>
<td headers="Frontal  q Value" class="gt_row gt_right">0.024</td>
<td headers="Frontal  Max fold change" class="gt_row gt_right">1.489</td>
<td headers="Frontal  Power" class="gt_row gt_right">0.568</td>
<td headers="Frontal  Highest mean condition" class="gt_row gt_left">2A</td>
<td headers="Frontal  Lowest mean condition" class="gt_row gt_left">1A</td>
<td headers="Frontal  highest_mean" class="gt_row gt_left">POS</td>
<td headers="Frontal  lowest_mean" class="gt_row gt_left">POS</td>
<td headers="Frontal  Mass" class="gt_row gt_right">101600</td>
<td headers="Frontal  Description" class="gt_row gt_left">Sodium/calcium exchanger 2 OS=Rattus norvegicus GN=Slc8a2 PE=1 SV=1</td></tr>
    <tr><td headers="Frontal  Accession" class="gt_row gt_left">TCPD_RAT;TCPD_MOUSE</td>
<td headers="Frontal  Peptide count" class="gt_row gt_right">12</td>
<td headers="Frontal  Unique peptides" class="gt_row gt_right">2</td>
<td headers="Frontal  Confidence score" class="gt_row gt_right">93.440</td>
<td headers="Frontal  Anova (p)" class="gt_row gt_right">0.049</td>
<td headers="Frontal  q Value" class="gt_row gt_right">0.024</td>
<td headers="Frontal  Max fold change" class="gt_row gt_right">2.945</td>
<td headers="Frontal  Power" class="gt_row gt_right">0.569</td>
<td headers="Frontal  Highest mean condition" class="gt_row gt_left">3A</td>
<td headers="Frontal  Lowest mean condition" class="gt_row gt_left">1A</td>
<td headers="Frontal  highest_mean" class="gt_row gt_left">POS</td>
<td headers="Frontal  lowest_mean" class="gt_row gt_left">POS</td>
<td headers="Frontal  Mass" class="gt_row gt_right">58576</td>
<td headers="Frontal  Description" class="gt_row gt_left">T-complex protein 1 subunit delta OS=Rattus norvegicus GN=Cct4 PE=1 SV=3</td></tr>
    <tr class="gt_group_heading_row">
      <th colspan="14" class="gt_group_heading" scope="colgroup" id="Hippocampus">Hippocampus</th>
    </tr>
    <tr class="gt_row_group_first"><td headers="Hippocampus  Accession" class="gt_row gt_left">KPCA_MOUSE;KPCA_RAT</td>
<td headers="Hippocampus  Peptide count" class="gt_row gt_right">10</td>
<td headers="Hippocampus  Unique peptides" class="gt_row gt_right">3</td>
<td headers="Hippocampus  Confidence score" class="gt_row gt_right">131.200</td>
<td headers="Hippocampus  Anova (p)" class="gt_row gt_right">0.003</td>
<td headers="Hippocampus  q Value" class="gt_row gt_right">1.000</td>
<td headers="Hippocampus  Max fold change" class="gt_row gt_right">3.684</td>
<td headers="Hippocampus  Power" class="gt_row gt_right">0.946</td>
<td headers="Hippocampus  Highest mean condition" class="gt_row gt_left">2B</td>
<td headers="Hippocampus  Lowest mean condition" class="gt_row gt_left">1B</td>
<td headers="Hippocampus  highest_mean" class="gt_row gt_left">NEG</td>
<td headers="Hippocampus  lowest_mean" class="gt_row gt_left">NEG</td>
<td headers="Hippocampus  Mass" class="gt_row gt_right">77943</td>
<td headers="Hippocampus  Description" class="gt_row gt_left">Protein kinase C alpha type OS=Mus musculus GN=Prkca PE=1 SV=3</td></tr>
    <tr><td headers="Hippocampus  Accession" class="gt_row gt_left">UBE2N_MOUSE</td>
<td headers="Hippocampus  Peptide count" class="gt_row gt_right">3</td>
<td headers="Hippocampus  Unique peptides" class="gt_row gt_right">1</td>
<td headers="Hippocampus  Confidence score" class="gt_row gt_right">110.320</td>
<td headers="Hippocampus  Anova (p)" class="gt_row gt_right">0.018</td>
<td headers="Hippocampus  q Value" class="gt_row gt_right">1.000</td>
<td headers="Hippocampus  Max fold change" class="gt_row gt_right">2.936</td>
<td headers="Hippocampus  Power" class="gt_row gt_right">0.763</td>
<td headers="Hippocampus  Highest mean condition" class="gt_row gt_left">2B</td>
<td headers="Hippocampus  Lowest mean condition" class="gt_row gt_left">1B</td>
<td headers="Hippocampus  highest_mean" class="gt_row gt_left">NEG</td>
<td headers="Hippocampus  lowest_mean" class="gt_row gt_left">NEG</td>
<td headers="Hippocampus  Mass" class="gt_row gt_right">17184</td>
<td headers="Hippocampus  Description" class="gt_row gt_left">Ubiquitin-conjugating enzyme E2 N OS=Mus musculus GN=Ube2n PE=1 SV=1</td></tr>
    <tr><td headers="Hippocampus  Accession" class="gt_row gt_left">CPNE3_MOUSE</td>
<td headers="Hippocampus  Peptide count" class="gt_row gt_right">2</td>
<td headers="Hippocampus  Unique peptides" class="gt_row gt_right">0</td>
<td headers="Hippocampus  Confidence score" class="gt_row gt_right">5.240</td>
<td headers="Hippocampus  Anova (p)" class="gt_row gt_right">0.026</td>
<td headers="Hippocampus  q Value" class="gt_row gt_right">1.000</td>
<td headers="Hippocampus  Max fold change" class="gt_row gt_right">2.709</td>
<td headers="Hippocampus  Power" class="gt_row gt_right">0.705</td>
<td headers="Hippocampus  Highest mean condition" class="gt_row gt_left">1B</td>
<td headers="Hippocampus  Lowest mean condition" class="gt_row gt_left">3B</td>
<td headers="Hippocampus  highest_mean" class="gt_row gt_left">NEG</td>
<td headers="Hippocampus  lowest_mean" class="gt_row gt_left">NEG</td>
<td headers="Hippocampus  Mass" class="gt_row gt_right">60288</td>
<td headers="Hippocampus  Description" class="gt_row gt_left">Copine-3 OS=Mus musculus GN=Cpne3 PE=1 SV=2</td></tr>
    <tr><td headers="Hippocampus  Accession" class="gt_row gt_left">CMC1_MOUSE</td>
<td headers="Hippocampus  Peptide count" class="gt_row gt_right">10</td>
<td headers="Hippocampus  Unique peptides" class="gt_row gt_right">3</td>
<td headers="Hippocampus  Confidence score" class="gt_row gt_right">383.050</td>
<td headers="Hippocampus  Anova (p)" class="gt_row gt_right">0.033</td>
<td headers="Hippocampus  q Value" class="gt_row gt_right">1.000</td>
<td headers="Hippocampus  Max fold change" class="gt_row gt_right">1.604</td>
<td headers="Hippocampus  Power" class="gt_row gt_right">0.665</td>
<td headers="Hippocampus  Highest mean condition" class="gt_row gt_left">3B</td>
<td headers="Hippocampus  Lowest mean condition" class="gt_row gt_left">1B</td>
<td headers="Hippocampus  highest_mean" class="gt_row gt_left">NEG</td>
<td headers="Hippocampus  lowest_mean" class="gt_row gt_left">NEG</td>
<td headers="Hippocampus  Mass" class="gt_row gt_right">74922</td>
<td headers="Hippocampus  Description" class="gt_row gt_left">Calcium-binding mitochondrial carrier protein Aralar1 OS=Mus musculus GN=Slc25a12 PE=1 SV=1</td></tr>
    <tr><td headers="Hippocampus  Accession" class="gt_row gt_left">ATPD_MOUSE</td>
<td headers="Hippocampus  Peptide count" class="gt_row gt_right">5</td>
<td headers="Hippocampus  Unique peptides" class="gt_row gt_right">2</td>
<td headers="Hippocampus  Confidence score" class="gt_row gt_right">173.480</td>
<td headers="Hippocampus  Anova (p)" class="gt_row gt_right">0.046</td>
<td headers="Hippocampus  q Value" class="gt_row gt_right">1.000</td>
<td headers="Hippocampus  Max fold change" class="gt_row gt_right">2.148</td>
<td headers="Hippocampus  Power" class="gt_row gt_right">0.600</td>
<td headers="Hippocampus  Highest mean condition" class="gt_row gt_left">2B</td>
<td headers="Hippocampus  Lowest mean condition" class="gt_row gt_left">3B</td>
<td headers="Hippocampus  highest_mean" class="gt_row gt_left">NEG</td>
<td headers="Hippocampus  lowest_mean" class="gt_row gt_left">NEG</td>
<td headers="Hippocampus  Mass" class="gt_row gt_right">17589</td>
<td headers="Hippocampus  Description" class="gt_row gt_left">ATP synthase subunit delta, mitochondrial OS=Mus musculus GN=Atp5d PE=1 SV=1</td></tr>
  </tbody>
  &#10;</table>
</div>

``` r
  # tab_style(
  #   style = cell_fill(color = "lightblue"),
  #   locations = cells_body(rows = log2FC < 0)
  # )
```

------------------------------------------------------------------------

# Part 2 - mapping & other stuff

``` r
# load in uniprot df
load(file = "input/mouse_ref_proteome_2025_12_13_54864x295.rda")
load(file = "input/rat_ref_proteome_2025_12_13_51872x295.rda")
```

``` r
# cols I want from uniprot

#colnames(mouse_ref_proteome_2025_12_13)

my_cols <- c("Entry", "Entry.Name", "Protein.names", "Gene.Names..primary.",
             "Gene.Names", "Gene.Names..synonym.",  "Annotation", "Reviewed",
             "Protein.existence", "Protein.families", "Length", "Organism",
             "Keywords", "Keyword.ID", "Function..CC.", "Miscellaneous..CC.", "Pathway",
             "GeneID", "Ensembl",
             "Gene.Ontology..GO.", "Gene.Ontology..cellular.component.",
             "Gene.Ontology..molecular.function.", "Gene.Ontology..biological.process.",
             "Gene.Ontology.IDs", "KEGG", "Reactome", "STRING",
             "Entry.version", "Date.of.creation", "Date.of.last.modification",
             "Date.of.last.sequence.modification"
             )
```

``` r
# MOUSE REF PROTEOME
mouse_ref <- mouse_ref_proteome_2025_12_13 %>% 
  dplyr::select(all_of(my_cols)) %>%
  rename(gene_names_primary = Gene.Names..primary.,
         gene_names_synonym = Gene.Names..synonym.) %>% 
  mutate(existence_ranking = case_when( # ranking level of evidence
    Protein.existence == "Evidence at protein level" ~ 5,
    Protein.existence == "Evidence at transcript level" ~ 4,
    Protein.existence == "Inferred from homology" ~ 3, 
    Protein.existence == "Predicted" ~ 2, 
    Protein.existence == "Uncertain" ~ 1
    ), .after = Protein.existence) %>%
  mutate(across(where(is.character), ~ na_if(.x, ""))) 


# need to run in separate line as mouse ref is input (otheriwise object not found )
mouse_ref <- mouse_ref %>% 
  mutate(occupied = rowSums(!is.na(mouse_ref)))

# 2025_12_15 - NOTE DID NOT GET MOUSE OTHER PROTEOME FROM UNIPROT BECAUSE IT SAYS CONTAIMINATED!
```

``` r
# RAT REF PROTEOME
rat_ref <- rat_ref_proteome_2025_12_13 %>% 
  dplyr::select(all_of(my_cols)) %>%
  rename(gene_names_primary = Gene.Names..primary.,
         gene_names_synonym = Gene.Names..synonym.) %>% 
  mutate(existence_ranking = case_when( # ranking level of evidence
    Protein.existence == "Evidence at protein level" ~ 5,
    Protein.existence == "Evidence at transcript level" ~ 4,
    Protein.existence == "Inferred from homology" ~ 3, 
    Protein.existence == "Predicted" ~ 2, 
    Protein.existence == "Uncertain" ~ 1
    ), .after = Protein.existence) %>%
  mutate(across(where(is.character), ~ na_if(.x, ""))) 

rat_ref <- rat_ref %>%
  mutate(occupied = rowSums(!is.na(rat_ref)))
```

``` r
combined_2 <- combined %>%
  # Prep
  mutate(row_global = row_number(), .before = Region) %>% # global row no.
  group_by(Region) %>% 
  mutate(row_region = row_number(), .after = Region) %>% # row no. by region
  
  # String split (some IDs are KPCA_MOUSE;KPCA_RAT)
  mutate(split_acc = str_split(Accession, ";"), .after = Accession) %>%
  tidyr::unnest(split_acc) %>%  # goes from 1954 to 2360
  ungroup()
```

``` r
rat_shaw_df_1 <- combined_2 %>%
  filter(str_detect(split_acc, "_RAT")) %>% # 692 proteins
  left_join(rat_ref, by = join_by(split_acc == Entry.Name)) %>%
  group_by(row_global) %>%
  filter(if(any(Reviewed == "reviewed", na.rm = TRUE)) {
    Reviewed == "reviewed"
  } else{TRUE}) %>% #692 
  slice_max(Annotation, n =1) %>%  #690
  slice_max(existence_ranking, n = 1) %>%  #690
  slice_max(occupied, n = 1, with_ties = FALSE) %>%  #688 | 686 no ties
  ungroup()

sum(is.na(rat_shaw_df_1$Entry)) # 12 unmapped (seem like non ref ones? or removed ones?)
```

    ## [1] 12

``` r
length(unique(rat_shaw_df_1$Entry)) #341
```

    ## [1] 338

``` r
length(unique(rat_shaw_df_1$row_global)) #686
```

    ## [1] 686

``` r
# rat_shaw_df_1 %>%
#   count(Entry) %>%
#   arrange(desc(n))
```

``` r
mouse_shaw_df_1 <- combined_2 %>%
  filter(str_detect(split_acc, "_MOUSE")) %>% #1390
  left_join(mouse_ref, by = join_by(split_acc == Entry.Name)) %>%
  rename(Entry.m = Entry, # rename so no doubleups when joining to rat
         Protein.names.m = Protein.names,
         gene_names_primary.m = gene_names_primary,
         Gene.Names.m = Gene.Names,
         gene_names_synonym.m = gene_names_synonym)

sum(is.na(mouse_shaw_df_1$Entry.m)) # 26 unmapped
```

    ## [1] 26

``` r
length(unique(mouse_shaw_df_1$Entry.m)) #683
```

    ## [1] 683

``` r
length(unique(mouse_shaw_df_1$row_global)) #1350
```

    ## [1] 1350

``` r
mouse_2_rat_shaw <- mouse_shaw_df_1 %>%  #1390 
  select(row_global:gene_names_synonym.m) %>% # peelback
  left_join(rat_ref, by = join_by(gene_names_primary.m == gene_names_primary),
            relationship = "many-to-many") %>%  #34448
  group_by(row_global) %>%
  filter(if(any(Reviewed == "reviewed", na.rm = TRUE)) {
    Reviewed == "reviewed"
  } else{TRUE}) %>%  #5188 
  slice_max(Annotation, n=1) %>%  #1828
  slice_max(existence_ranking, n = 1) %>%  #1720
  slice_max(occupied, n=1, with_ties = FALSE ) %>%  #1402 | no ties 1350 
  ungroup()
```

### Make combined mapped df

Ignore

``` r
combined_map_df <- bind_rows(list(rat_entry_2_rat_acc = rat_shaw_df_1,
                                  mouse_entry_2_rat_gene = mouse_2_rat_shaw),
                             .id = "map_method_df") %>% #2036
  group_by(row_global) %>%
  filter(if(any(map_method_df == "rat_entry_2_rat_acc", na.rm = TRUE)) {
    map_method_df == "rat_entry_2_rat_acc"
  } else{TRUE}) %>% # 1792 
  ungroup() %>%
  mutate(dir = ifelse(`Max fold change` > 0, "Down", "Up"),
         .after = `Max fold change`) %>% # note opposite btw! since alll were up in NEG
  rename(Region_ = Region)
```

### 210 entries from other rodents!

``` r
# all other rodents! - so many different rodents!
# CRIGR_ chinese hamster..., CAVPO_ guinea pig...,
other_unmapped_species <- combined_2 %>%
  filter(!str_detect(Accession, "_RAT|_MOUSE")) %>% # 210 rows with n
  arrange(`Anova (p)`)

head(other_unmapped_species) # only ODP2_MESAU - part of pyruvate dehydrogenase complex
```

    ## # A tibble: 6 × 18
    ##   row_global Region      row_region Accession          split_acc `Peptide count`
    ##        <int> <chr>            <int> <chr>              <chr>     <chr>          
    ## 1          1 Frontal              1 ODP2_MESAU         ODP2_MES… 7              
    ## 2         50 Frontal             50 COX2_GERGE         COX2_GER… 7              
    ## 3         51 Frontal             51 GNAI2_CAVPO        GNAI2_CA… 7              
    ## 4         63 Frontal             63 COX2_RHAPU;COX2_R… COX2_RHA… 7              
    ## 5         63 Frontal             63 COX2_RHAPU;COX2_R… COX2_RHY… 7              
    ## 6        986 Hippocampus          9 GNAI2_CAVPO        GNAI2_CA… 7              
    ## # ℹ 12 more variables: `Unique peptides` <chr>, `Confidence score` <dbl>,
    ## #   `Anova (p)` <dbl>, `q Value` <dbl>, `Max fold change` <dbl>, Power <dbl>,
    ## #   `Highest mean condition` <chr>, `Lowest mean condition` <chr>,
    ## #   highest_mean <chr>, lowest_mean <chr>, Mass <chr>, Description <chr>

------------------------------------------------------------------------

## Lit Review Table

- <https://hbctraining.github.io/Intro-to-R/lessons/basic_plots_in_r.html>
- <https://rstudio-pubs-static.s3.amazonaws.com/7953_4e3efd5b9415444ca065b1167862c349.html>

``` r
colnames(combined_map_df)
```

    ##  [1] "map_method_df"                      "row_global"                        
    ##  [3] "Region_"                            "row_region"                        
    ##  [5] "Accession"                          "split_acc"                         
    ##  [7] "Peptide count"                      "Unique peptides"                   
    ##  [9] "Confidence score"                   "Anova (p)"                         
    ## [11] "q Value"                            "Max fold change"                   
    ## [13] "dir"                                "Power"                             
    ## [15] "Highest mean condition"             "Lowest mean condition"             
    ## [17] "highest_mean"                       "lowest_mean"                       
    ## [19] "Mass"                               "Description"                       
    ## [21] "Entry"                              "Protein.names"                     
    ## [23] "gene_names_primary"                 "Gene.Names"                        
    ## [25] "gene_names_synonym"                 "Annotation"                        
    ## [27] "Reviewed"                           "Protein.existence"                 
    ## [29] "existence_ranking"                  "Protein.families"                  
    ## [31] "Length"                             "Organism"                          
    ## [33] "Keywords"                           "Keyword.ID"                        
    ## [35] "Function..CC."                      "Miscellaneous..CC."                
    ## [37] "Pathway"                            "GeneID"                            
    ## [39] "Ensembl"                            "Gene.Ontology..GO."                
    ## [41] "Gene.Ontology..cellular.component." "Gene.Ontology..molecular.function."
    ## [43] "Gene.Ontology..biological.process." "Gene.Ontology.IDs"                 
    ## [45] "KEGG"                               "Reactome"                          
    ## [47] "STRING"                             "Entry.version"                     
    ## [49] "Date.of.creation"                   "Date.of.last.modification"         
    ## [51] "Date.of.last.sequence.modification" "occupied"                          
    ## [53] "Entry.m"                            "Protein.names.m"                   
    ## [55] "gene_names_primary.m"               "Gene.Names.m"                      
    ## [57] "gene_names_synonym.m"               "Entry.Name"

``` r
plot(combined_map_df$`Anova (p)`) # not super extra most 0.5 or -0.5
```

![](Shaw_2020_mice_proteomics_POSvsNEG_data_files/figure-gfm/unnamed-chunk-27-1.png)<!-- -->

``` r
plot(combined_map_df$`q Value`)
```

![](Shaw_2020_mice_proteomics_POSvsNEG_data_files/figure-gfm/unnamed-chunk-27-2.png)<!-- -->

``` r
plot(combined_map_df$Power)
```

![](Shaw_2020_mice_proteomics_POSvsNEG_data_files/figure-gfm/unnamed-chunk-27-3.png)<!-- -->

``` r
plot(combined_map_df$`Max fold change`)
```

![](Shaw_2020_mice_proteomics_POSvsNEG_data_files/figure-gfm/unnamed-chunk-27-4.png)<!-- -->

``` r
plot(combined_map_df$`Max fold change`, ylim=c(0,7)) # zoomed in
```

![](Shaw_2020_mice_proteomics_POSvsNEG_data_files/figure-gfm/unnamed-chunk-27-5.png)<!-- -->

``` r
# # no colors???
# # plot(`Anova (p)` ~ `Max fold change`, data = combined_map_df,
# #      # # main = "", xlab = "", ylab = "",
# #      col = c("pink","orange")[combined_map_df$Region], pch=16, cex=2.0)
# 
# plot(combined_map_df$`Anova (p)`, combined_map_df$`Max fold change`,
#      col = combined_map_df$Region) 
```

``` r
# for lit review table 
lit_review_proteomics <- combined_map_df %>%
  mutate(
    space = " ",
    Year = 2020,
    Author = 'Shaw',
    Animal = 'Mice',
    Condition = ' ',
    Region = Region_ , 
    Whole = 'Yes',
    Exp_type = 'Protein',
    Method = 'PROTEOMICS_LCMS_LFQ_DDA_OS',
    Used_name = paste0(Accession, ' - (', Description,')'),
    Mod = " ",
    Mod_type = " ",
    Gene_name = ifelse(map_method_df == "rat_entry_2_rat_acc",
                       gene_names_primary, # if rat 2 rat 
                       gene_names_primary.m), # if mouse 2 rat (other wise gives NA)
    Comparison = 4, # POS VS NEG
    In_EE = dir,
    Fold_change = paste0("-", round(`Max fold change`, 2)), # prob same as log2FC?
    Arrows = case_when( # dplyr case when, instead of ifelse
      `Max fold change` > 2.5 ~ '+++',
      `Max fold change` > 1 ~ '++',
      `Max fold change` > .05 ~ '+',
      `Max fold change` > -.05 ~ '=',
      `Max fold change` > -1 ~ '-',
      `Max fold change` > -2.5 ~ '--',
      `Max fold change` <= -2.5 ~ '---',
      ),
    Significant = ifelse(`Anova (p)` <= 0.05, 'yes', 'no'),
    pvalue = `Anova (p)`, # not they have q val (FDR which was used) & etc, 
    N_tested = '5',
    Stats_adj = 'unadjusted',
    Note = '',
    Full_name = Protein.names,
    Uniprot_rat_accession = Entry,
    # Uniprot_entry_name = Entry.Name,
    Uniprot_entry_name = ifelse(map_method_df == "rat_entry_2_rat_acc",
                       split_acc, # if rat 2 rat (other NA) since was used for joining
                       Entry.Name), # else Entry.Name (normally here)
    Uniprot_link = ifelse(!is.na(Entry),
                          paste0('https://www.uniprot.org/uniprotkb/', Entry,'/entry'),
                          ' '),
    Wikipedia_link = ifelse(!is.na(Entry),
                            paste0('https://en.wikipedia.org/wiki/',
                                   `gene_names_primary`), ' '),
    
    Curation_method = case_when(
      Reviewed == "reviewed" & map_method_df == "rat_entry_2_rat_acc"
        ~ "auto - ref & reviewed (rat_entry_2_rat_acc)",
      
      Reviewed == "unreviewed" & map_method_df == "rat_entry_2_rat_acc"
        ~ "auto - ref & unreviewed (rat_entry_2_rat_acc)",
      
      Reviewed == "reviewed" & map_method_df == "mouse_entry_2_rat_gene"
        ~ "auto - ref & reviewed (mouse_entry_2_rat_gene_primary)",
      
      Reviewed == "unreviewed" & map_method_df == "mouse_entry_2_rat_gene"
        ~ "auto - ref proteome (mouse_entry_2_rat_gene_primary)",
      
      is.na(Entry.Name)
        ~ "manual - na"
      ),
    Class_group = ' ',
    Other_link = ' ',
    Other_comments = 'misc is combined EntryNameAcc_Description from the provided df',
    misc = paste0(Accession, '_', Description)
  ) %>%
  arrange(Region, desc(In_EE), Gene_name)
```

``` r
# # write out | #1792
# write_excel_csv(lit_review_proteomics,
#     "output/2025_12_29_Shaw_2020_mice_hip_fc_proteomics_lit_review_allratmouse_v1.csv")
```

#### Significant Only

``` r
lit_review_proteomics_sig <- lit_review_proteomics %>%
  filter(`Anova (p)` <= 0.05) # 49
```

``` r
# #write out
# write_excel_csv(lit_review_proteomics_sig,
#     "output/2025_12_29_Shaw_2020_mice_hip_fc_proteomics_lit_review_sig_ramo_v1.csv")
```

## gt table lit mapped

``` r
lit_review_proteomics_sig |> 
  
  dplyr::select(c(Region, In_EE, Gene_name, Entry, Uniprot_entry_name, Full_name, 
                  Accession, Description,
                  Fold_change, `Anova (p)`, Reviewed, Protein.families)) |>
  
  gt(
    groupname_col = "Region",
    rowname_col = "In_EE"
  ) |>
  
  cols_label(
    Gene_name = "Entry",
    Uniprot_entry_name = "Name",
    Full_name = "Protein Description",
    # misc = "Arroyo's Identifiers"
    Accession = "Shaw Acc",
    Description = "Shaw Desc",
  ) |>
  
  fmt_number(decimals = 2) |>
  
  tab_header(
    title = "Shaw et al., 2020 (Mice Hippocampus & Frontal Pole/Lobe - POS(EE) vs. NEG)",
    subtitle = "44 DEG (Front), 5 DEG (Hip) - all down in POS(EE) p<=.05, data is a bit odd/confusing (removed excess rodents) + some missing data? look at note at the start of the knitted", 
  ) |>
  
  tab_style( 
    style = list(
      # cell_fill(color = "#F9E3D6"),
      cell_text(style = "italic")
      ),
    locations = cells_body(
      columns = Gene_name,
      # rows = currency < 100
    )
  ) |>
  
  tab_style( 
    style = list(
      # cell_fill(color = "#F9E3D6"),
      cell_text(size = 8)
      ),
    locations = cells_body(
      columns = c(Accession, Description, Reviewed)
      # rows = currency < 100
    )
  ) |>
  
  tab_options(heading.title.font.size = 16.5,
              table.font.size=11.5,
              # stub.font.weight = "bold",
              # stub.font.size = 13,
              row_group.font.weight = "bold",
              row_group.font.size = 14,
              column_labels.font.size = 12.5,
              column_labels.font.weight = "bold") #|>
```

<div id="xjlundoxcj" style="padding-left:0px;padding-right:0px;padding-top:10px;padding-bottom:10px;overflow-x:auto;overflow-y:auto;width:auto;height:auto;">
<style>#xjlundoxcj table {
  font-family: system-ui, 'Segoe UI', Roboto, Helvetica, Arial, sans-serif, 'Apple Color Emoji', 'Segoe UI Emoji', 'Segoe UI Symbol', 'Noto Color Emoji';
  -webkit-font-smoothing: antialiased;
  -moz-osx-font-smoothing: grayscale;
}
&#10;#xjlundoxcj thead, #xjlundoxcj tbody, #xjlundoxcj tfoot, #xjlundoxcj tr, #xjlundoxcj td, #xjlundoxcj th {
  border-style: none;
}
&#10;#xjlundoxcj p {
  margin: 0;
  padding: 0;
}
&#10;#xjlundoxcj .gt_table {
  display: table;
  border-collapse: collapse;
  line-height: normal;
  margin-left: auto;
  margin-right: auto;
  color: #333333;
  font-size: 11.5px;
  font-weight: normal;
  font-style: normal;
  background-color: #FFFFFF;
  width: auto;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #A8A8A8;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #A8A8A8;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
}
&#10;#xjlundoxcj .gt_caption {
  padding-top: 4px;
  padding-bottom: 4px;
}
&#10;#xjlundoxcj .gt_title {
  color: #333333;
  font-size: 16.5px;
  font-weight: initial;
  padding-top: 4px;
  padding-bottom: 4px;
  padding-left: 5px;
  padding-right: 5px;
  border-bottom-color: #FFFFFF;
  border-bottom-width: 0;
}
&#10;#xjlundoxcj .gt_subtitle {
  color: #333333;
  font-size: 85%;
  font-weight: initial;
  padding-top: 3px;
  padding-bottom: 5px;
  padding-left: 5px;
  padding-right: 5px;
  border-top-color: #FFFFFF;
  border-top-width: 0;
}
&#10;#xjlundoxcj .gt_heading {
  background-color: #FFFFFF;
  text-align: center;
  border-bottom-color: #FFFFFF;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
}
&#10;#xjlundoxcj .gt_bottom_border {
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}
&#10;#xjlundoxcj .gt_col_headings {
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
}
&#10;#xjlundoxcj .gt_col_heading {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 12.5px;
  font-weight: bold;
  text-transform: inherit;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: bottom;
  padding-top: 5px;
  padding-bottom: 6px;
  padding-left: 5px;
  padding-right: 5px;
  overflow-x: hidden;
}
&#10;#xjlundoxcj .gt_column_spanner_outer {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 12.5px;
  font-weight: bold;
  text-transform: inherit;
  padding-top: 0;
  padding-bottom: 0;
  padding-left: 4px;
  padding-right: 4px;
}
&#10;#xjlundoxcj .gt_column_spanner_outer:first-child {
  padding-left: 0;
}
&#10;#xjlundoxcj .gt_column_spanner_outer:last-child {
  padding-right: 0;
}
&#10;#xjlundoxcj .gt_column_spanner {
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  vertical-align: bottom;
  padding-top: 5px;
  padding-bottom: 5px;
  overflow-x: hidden;
  display: inline-block;
  width: 100%;
}
&#10;#xjlundoxcj .gt_spanner_row {
  border-bottom-style: hidden;
}
&#10;#xjlundoxcj .gt_group_heading {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  color: #333333;
  background-color: #FFFFFF;
  font-size: 14px;
  font-weight: bold;
  text-transform: inherit;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: middle;
  text-align: left;
}
&#10;#xjlundoxcj .gt_empty_group_heading {
  padding: 0.5px;
  color: #333333;
  background-color: #FFFFFF;
  font-size: 14px;
  font-weight: bold;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  vertical-align: middle;
}
&#10;#xjlundoxcj .gt_from_md > :first-child {
  margin-top: 0;
}
&#10;#xjlundoxcj .gt_from_md > :last-child {
  margin-bottom: 0;
}
&#10;#xjlundoxcj .gt_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  margin: 10px;
  border-top-style: solid;
  border-top-width: 1px;
  border-top-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: middle;
  overflow-x: hidden;
}
&#10;#xjlundoxcj .gt_stub {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  text-transform: inherit;
  border-right-style: solid;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
  padding-left: 5px;
  padding-right: 5px;
}
&#10;#xjlundoxcj .gt_stub_row_group {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  text-transform: inherit;
  border-right-style: solid;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
  padding-left: 5px;
  padding-right: 5px;
  vertical-align: top;
}
&#10;#xjlundoxcj .gt_row_group_first td {
  border-top-width: 2px;
}
&#10;#xjlundoxcj .gt_row_group_first th {
  border-top-width: 2px;
}
&#10;#xjlundoxcj .gt_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}
&#10;#xjlundoxcj .gt_first_summary_row {
  border-top-style: solid;
  border-top-color: #D3D3D3;
}
&#10;#xjlundoxcj .gt_first_summary_row.thick {
  border-top-width: 2px;
}
&#10;#xjlundoxcj .gt_last_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}
&#10;#xjlundoxcj .gt_grand_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}
&#10;#xjlundoxcj .gt_first_grand_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-top-style: double;
  border-top-width: 6px;
  border-top-color: #D3D3D3;
}
&#10;#xjlundoxcj .gt_last_grand_summary_row_top {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-bottom-style: double;
  border-bottom-width: 6px;
  border-bottom-color: #D3D3D3;
}
&#10;#xjlundoxcj .gt_striped {
  background-color: rgba(128, 128, 128, 0.05);
}
&#10;#xjlundoxcj .gt_table_body {
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}
&#10;#xjlundoxcj .gt_footnotes {
  color: #333333;
  background-color: #FFFFFF;
  border-bottom-style: none;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
}
&#10;#xjlundoxcj .gt_footnote {
  margin: 0px;
  font-size: 90%;
  padding-top: 4px;
  padding-bottom: 4px;
  padding-left: 5px;
  padding-right: 5px;
}
&#10;#xjlundoxcj .gt_sourcenotes {
  color: #333333;
  background-color: #FFFFFF;
  border-bottom-style: none;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
}
&#10;#xjlundoxcj .gt_sourcenote {
  font-size: 90%;
  padding-top: 4px;
  padding-bottom: 4px;
  padding-left: 5px;
  padding-right: 5px;
}
&#10;#xjlundoxcj .gt_left {
  text-align: left;
}
&#10;#xjlundoxcj .gt_center {
  text-align: center;
}
&#10;#xjlundoxcj .gt_right {
  text-align: right;
  font-variant-numeric: tabular-nums;
}
&#10;#xjlundoxcj .gt_font_normal {
  font-weight: normal;
}
&#10;#xjlundoxcj .gt_font_bold {
  font-weight: bold;
}
&#10;#xjlundoxcj .gt_font_italic {
  font-style: italic;
}
&#10;#xjlundoxcj .gt_super {
  font-size: 65%;
}
&#10;#xjlundoxcj .gt_footnote_marks {
  font-size: 75%;
  vertical-align: 0.4em;
  position: initial;
}
&#10;#xjlundoxcj .gt_asterisk {
  font-size: 100%;
  vertical-align: 0;
}
&#10;#xjlundoxcj .gt_indent_1 {
  text-indent: 5px;
}
&#10;#xjlundoxcj .gt_indent_2 {
  text-indent: 10px;
}
&#10;#xjlundoxcj .gt_indent_3 {
  text-indent: 15px;
}
&#10;#xjlundoxcj .gt_indent_4 {
  text-indent: 20px;
}
&#10;#xjlundoxcj .gt_indent_5 {
  text-indent: 25px;
}
&#10;#xjlundoxcj .katex-display {
  display: inline-flex !important;
  margin-bottom: 0.75em !important;
}
&#10;#xjlundoxcj div.Reactable > div.rt-table > div.rt-thead > div.rt-tr.rt-tr-group-header > div.rt-th-group:after {
  height: 0px !important;
}
</style>
<table class="gt_table" data-quarto-disable-processing="false" data-quarto-bootstrap="false">
  <thead>
    <tr class="gt_heading">
      <td colspan="11" class="gt_heading gt_title gt_font_normal" style>Shaw et al., 2020 (Mice Hippocampus &amp; Frontal Pole/Lobe - POS(EE) vs. NEG)</td>
    </tr>
    <tr class="gt_heading">
      <td colspan="11" class="gt_heading gt_subtitle gt_font_normal gt_bottom_border" style>44 DEG (Front), 5 DEG (Hip) - all down in POS(EE) p&lt;=.05, data is a bit odd/confusing (removed excess rodents) + some missing data? look at note at the start of the knitted</td>
    </tr>
    <tr class="gt_col_headings">
      <th class="gt_col_heading gt_columns_bottom_border gt_left" rowspan="1" colspan="1" scope="col" id="a::stub"></th>
      <th class="gt_col_heading gt_columns_bottom_border gt_left" rowspan="1" colspan="1" scope="col" id="Gene_name">Entry</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_left" rowspan="1" colspan="1" scope="col" id="Entry">Entry</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_left" rowspan="1" colspan="1" scope="col" id="Uniprot_entry_name">Name</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_left" rowspan="1" colspan="1" scope="col" id="Full_name">Protein Description</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_left" rowspan="1" colspan="1" scope="col" id="Accession">Shaw Acc</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_left" rowspan="1" colspan="1" scope="col" id="Description">Shaw Desc</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="Fold_change">Fold_change</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="Anova-(p)">Anova (p)</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_left" rowspan="1" colspan="1" scope="col" id="Reviewed">Reviewed</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_left" rowspan="1" colspan="1" scope="col" id="Protein.families">Protein.families</th>
    </tr>
  </thead>
  <tbody class="gt_table_body">
    <tr class="gt_group_heading_row">
      <th colspan="11" class="gt_group_heading" scope="colgroup" id="Frontal">Frontal</th>
    </tr>
    <tr class="gt_row_group_first"><th id="stub_1_1" scope="row" class="gt_row gt_left gt_stub">Down</th>
<td headers="Frontal stub_1_1 Gene_name" class="gt_row gt_left" style="font-style: italic;">Acat1</td>
<td headers="Frontal stub_1_1 Entry" class="gt_row gt_left">P17764</td>
<td headers="Frontal stub_1_1 Uniprot_entry_name" class="gt_row gt_left">THIL_RAT</td>
<td headers="Frontal stub_1_1 Full_name" class="gt_row gt_left">Acetyl-CoA acetyltransferase, mitochondrial (EC 2.3.1.9) (Acetoacetyl-CoA thiolase)</td>
<td headers="Frontal stub_1_1 Accession" class="gt_row gt_left" style="font-size: 8;">THIL_RAT</td>
<td headers="Frontal stub_1_1 Description" class="gt_row gt_left" style="font-size: 8;">Acetyl-CoA acetyltransferase, mitochondrial OS=Rattus norvegicus GN=Acat1 PE=1 SV=1</td>
<td headers="Frontal stub_1_1 Fold_change" class="gt_row gt_right">-1.92</td>
<td headers="Frontal stub_1_1 Anova (p)" class="gt_row gt_right">0.01</td>
<td headers="Frontal stub_1_1 Reviewed" class="gt_row gt_left" style="font-size: 8;">reviewed</td>
<td headers="Frontal stub_1_1 Protein.families" class="gt_row gt_left">Thiolase-like superfamily, Thiolase family</td></tr>
    <tr><th id="stub_1_2" scope="row" class="gt_row gt_left gt_stub">Down</th>
<td headers="Frontal stub_1_2 Gene_name" class="gt_row gt_left" style="font-style: italic;">Acot7</td>
<td headers="Frontal stub_1_2 Entry" class="gt_row gt_left">Q64559</td>
<td headers="Frontal stub_1_2 Uniprot_entry_name" class="gt_row gt_left">BACH_RAT</td>
<td headers="Frontal stub_1_2 Full_name" class="gt_row gt_left">Cytosolic acyl coenzyme A thioester hydrolase (EC 3.1.2.2) (ACH1) (ACT) (Acyl-CoA thioesterase 7) (Brain acyl-CoA hydrolase) (BACH) (CTE-IIa) (CTE-IIb) (CTE-II) (LACH1) (Long chain acyl-CoA thioester hydrolase) (MTE-II)</td>
<td headers="Frontal stub_1_2 Accession" class="gt_row gt_left" style="font-size: 8;">BACH_RAT</td>
<td headers="Frontal stub_1_2 Description" class="gt_row gt_left" style="font-size: 8;">Cytosolic acyl coenzyme A thioester hydrolase OS=Rattus norvegicus GN=Acot7 PE=1 SV=4</td>
<td headers="Frontal stub_1_2 Fold_change" class="gt_row gt_right">-1.71</td>
<td headers="Frontal stub_1_2 Anova (p)" class="gt_row gt_right">0.05</td>
<td headers="Frontal stub_1_2 Reviewed" class="gt_row gt_left" style="font-size: 8;">reviewed</td>
<td headers="Frontal stub_1_2 Protein.families" class="gt_row gt_left">NA</td></tr>
    <tr><th id="stub_1_3" scope="row" class="gt_row gt_left gt_stub">Down</th>
<td headers="Frontal stub_1_3 Gene_name" class="gt_row gt_left" style="font-style: italic;">Acot7</td>
<td headers="Frontal stub_1_3 Entry" class="gt_row gt_left">Q64559</td>
<td headers="Frontal stub_1_3 Uniprot_entry_name" class="gt_row gt_left">BACH_RAT</td>
<td headers="Frontal stub_1_3 Full_name" class="gt_row gt_left">Cytosolic acyl coenzyme A thioester hydrolase (EC 3.1.2.2) (ACH1) (ACT) (Acyl-CoA thioesterase 7) (Brain acyl-CoA hydrolase) (BACH) (CTE-IIa) (CTE-IIb) (CTE-II) (LACH1) (Long chain acyl-CoA thioester hydrolase) (MTE-II)</td>
<td headers="Frontal stub_1_3 Accession" class="gt_row gt_left" style="font-size: 8;">BACH_MOUSE</td>
<td headers="Frontal stub_1_3 Description" class="gt_row gt_left" style="font-size: 8;">Cytosolic acyl coenzyme A thioester hydrolase OS=Mus musculus GN=Acot7 PE=1 SV=2</td>
<td headers="Frontal stub_1_3 Fold_change" class="gt_row gt_right">-4.02</td>
<td headers="Frontal stub_1_3 Anova (p)" class="gt_row gt_right">0.04</td>
<td headers="Frontal stub_1_3 Reviewed" class="gt_row gt_left" style="font-size: 8;">reviewed</td>
<td headers="Frontal stub_1_3 Protein.families" class="gt_row gt_left">NA</td></tr>
    <tr><th id="stub_1_4" scope="row" class="gt_row gt_left gt_stub">Down</th>
<td headers="Frontal stub_1_4 Gene_name" class="gt_row gt_left" style="font-style: italic;">Actn3</td>
<td headers="Frontal stub_1_4 Entry" class="gt_row gt_left">Q8R4I6</td>
<td headers="Frontal stub_1_4 Uniprot_entry_name" class="gt_row gt_left">Q8R4I6_RAT</td>
<td headers="Frontal stub_1_4 Full_name" class="gt_row gt_left">Actinin alpha 3</td>
<td headers="Frontal stub_1_4 Accession" class="gt_row gt_left" style="font-size: 8;">ACTN3_MOUSE</td>
<td headers="Frontal stub_1_4 Description" class="gt_row gt_left" style="font-size: 8;">Alpha-actinin-3 OS=Mus musculus GN=Actn3 PE=2 SV=1</td>
<td headers="Frontal stub_1_4 Fold_change" class="gt_row gt_right">-4.1</td>
<td headers="Frontal stub_1_4 Anova (p)" class="gt_row gt_right">0.04</td>
<td headers="Frontal stub_1_4 Reviewed" class="gt_row gt_left" style="font-size: 8;">unreviewed</td>
<td headers="Frontal stub_1_4 Protein.families" class="gt_row gt_left">Alpha-actinin family</td></tr>
    <tr><th id="stub_1_5" scope="row" class="gt_row gt_left gt_stub">Down</th>
<td headers="Frontal stub_1_5 Gene_name" class="gt_row gt_left" style="font-style: italic;">Add2</td>
<td headers="Frontal stub_1_5 Entry" class="gt_row gt_left">Q05764</td>
<td headers="Frontal stub_1_5 Uniprot_entry_name" class="gt_row gt_left">ADDB_RAT</td>
<td headers="Frontal stub_1_5 Full_name" class="gt_row gt_left">Beta-adducin (Adducin-63) (Erythrocyte adducin subunit beta)</td>
<td headers="Frontal stub_1_5 Accession" class="gt_row gt_left" style="font-size: 8;">ADDB_MOUSE</td>
<td headers="Frontal stub_1_5 Description" class="gt_row gt_left" style="font-size: 8;">Beta-adducin OS=Mus musculus GN=Add2 PE=1 SV=4</td>
<td headers="Frontal stub_1_5 Fold_change" class="gt_row gt_right">-1.77</td>
<td headers="Frontal stub_1_5 Anova (p)" class="gt_row gt_right">0.02</td>
<td headers="Frontal stub_1_5 Reviewed" class="gt_row gt_left" style="font-size: 8;">reviewed</td>
<td headers="Frontal stub_1_5 Protein.families" class="gt_row gt_left">Aldolase class II family, Adducin subfamily</td></tr>
    <tr><th id="stub_1_6" scope="row" class="gt_row gt_left gt_stub">Down</th>
<td headers="Frontal stub_1_6 Gene_name" class="gt_row gt_left" style="font-style: italic;">Aprt</td>
<td headers="Frontal stub_1_6 Entry" class="gt_row gt_left">P36972</td>
<td headers="Frontal stub_1_6 Uniprot_entry_name" class="gt_row gt_left">APT_RAT</td>
<td headers="Frontal stub_1_6 Full_name" class="gt_row gt_left">Adenine phosphoribosyltransferase (APRT) (EC 2.4.2.7)</td>
<td headers="Frontal stub_1_6 Accession" class="gt_row gt_left" style="font-size: 8;">APT_MOUSE</td>
<td headers="Frontal stub_1_6 Description" class="gt_row gt_left" style="font-size: 8;">Adenine phosphoribosyltransferase OS=Mus musculus GN=Aprt PE=1 SV=2</td>
<td headers="Frontal stub_1_6 Fold_change" class="gt_row gt_right">-3.51</td>
<td headers="Frontal stub_1_6 Anova (p)" class="gt_row gt_right">0.00</td>
<td headers="Frontal stub_1_6 Reviewed" class="gt_row gt_left" style="font-size: 8;">reviewed</td>
<td headers="Frontal stub_1_6 Protein.families" class="gt_row gt_left">Purine/pyrimidine phosphoribosyltransferase family</td></tr>
    <tr><th id="stub_1_7" scope="row" class="gt_row gt_left gt_stub">Down</th>
<td headers="Frontal stub_1_7 Gene_name" class="gt_row gt_left" style="font-style: italic;">Aqp4</td>
<td headers="Frontal stub_1_7 Entry" class="gt_row gt_left">P47863</td>
<td headers="Frontal stub_1_7 Uniprot_entry_name" class="gt_row gt_left">AQP4_RAT</td>
<td headers="Frontal stub_1_7 Full_name" class="gt_row gt_left">Aquaporin-4 (AQP-4) (Mercurial-insensitive water channel) (MIWC) (WCH4)</td>
<td headers="Frontal stub_1_7 Accession" class="gt_row gt_left" style="font-size: 8;">AQP4_RAT;AQP4_MOUSE</td>
<td headers="Frontal stub_1_7 Description" class="gt_row gt_left" style="font-size: 8;">Aquaporin-4 OS=Rattus norvegicus GN=Aqp4 PE=1 SV=1</td>
<td headers="Frontal stub_1_7 Fold_change" class="gt_row gt_right">-2.1</td>
<td headers="Frontal stub_1_7 Anova (p)" class="gt_row gt_right">0.02</td>
<td headers="Frontal stub_1_7 Reviewed" class="gt_row gt_left" style="font-size: 8;">reviewed</td>
<td headers="Frontal stub_1_7 Protein.families" class="gt_row gt_left">MIP/aquaporin (TC 1.A.8) family</td></tr>
    <tr><th id="stub_1_8" scope="row" class="gt_row gt_left gt_stub">Down</th>
<td headers="Frontal stub_1_8 Gene_name" class="gt_row gt_left" style="font-style: italic;">Atp6v1d</td>
<td headers="Frontal stub_1_8 Entry" class="gt_row gt_left">Q6P503</td>
<td headers="Frontal stub_1_8 Uniprot_entry_name" class="gt_row gt_left">Q6P503_RAT</td>
<td headers="Frontal stub_1_8 Full_name" class="gt_row gt_left">V-type proton ATPase subunit D (V-ATPase 28 kDa accessory protein) (V-type proton ATPase subunit d) (Vacuolar proton pump subunit D)</td>
<td headers="Frontal stub_1_8 Accession" class="gt_row gt_left" style="font-size: 8;">VATD_MOUSE</td>
<td headers="Frontal stub_1_8 Description" class="gt_row gt_left" style="font-size: 8;">V-type proton ATPase subunit D OS=Mus musculus GN=Atp6v1d PE=1 SV=1</td>
<td headers="Frontal stub_1_8 Fold_change" class="gt_row gt_right">-2.37</td>
<td headers="Frontal stub_1_8 Anova (p)" class="gt_row gt_right">0.03</td>
<td headers="Frontal stub_1_8 Reviewed" class="gt_row gt_left" style="font-size: 8;">unreviewed</td>
<td headers="Frontal stub_1_8 Protein.families" class="gt_row gt_left">V-ATPase D subunit family</td></tr>
    <tr><th id="stub_1_9" scope="row" class="gt_row gt_left gt_stub">Down</th>
<td headers="Frontal stub_1_9 Gene_name" class="gt_row gt_left" style="font-style: italic;">Atp6v1f</td>
<td headers="Frontal stub_1_9 Entry" class="gt_row gt_left">P50408</td>
<td headers="Frontal stub_1_9 Uniprot_entry_name" class="gt_row gt_left">VATF_RAT</td>
<td headers="Frontal stub_1_9 Full_name" class="gt_row gt_left">V-type proton ATPase subunit F (V-ATPase subunit F) (V-ATPase 14 kDa subunit) (Vacuolar proton pump subunit F)</td>
<td headers="Frontal stub_1_9 Accession" class="gt_row gt_left" style="font-size: 8;">VATF_MOUSE</td>
<td headers="Frontal stub_1_9 Description" class="gt_row gt_left" style="font-size: 8;">V-type proton ATPase subunit F OS=Mus musculus GN=Atp6v1f PE=1 SV=2</td>
<td headers="Frontal stub_1_9 Fold_change" class="gt_row gt_right">-1.88</td>
<td headers="Frontal stub_1_9 Anova (p)" class="gt_row gt_right">0.04</td>
<td headers="Frontal stub_1_9 Reviewed" class="gt_row gt_left" style="font-size: 8;">reviewed</td>
<td headers="Frontal stub_1_9 Protein.families" class="gt_row gt_left">V-ATPase F subunit family</td></tr>
    <tr><th id="stub_1_10" scope="row" class="gt_row gt_left gt_stub">Down</th>
<td headers="Frontal stub_1_10 Gene_name" class="gt_row gt_left" style="font-style: italic;">Atp6v1h</td>
<td headers="Frontal stub_1_10 Entry" class="gt_row gt_left">A0A8I5ZWI2</td>
<td headers="Frontal stub_1_10 Uniprot_entry_name" class="gt_row gt_left">A0A8I5ZWI2_RAT</td>
<td headers="Frontal stub_1_10 Full_name" class="gt_row gt_left">V-type proton ATPase subunit H</td>
<td headers="Frontal stub_1_10 Accession" class="gt_row gt_left" style="font-size: 8;">VATH_MOUSE</td>
<td headers="Frontal stub_1_10 Description" class="gt_row gt_left" style="font-size: 8;">V-type proton ATPase subunit H OS=Mus musculus GN=Atp6v1h PE=1 SV=1</td>
<td headers="Frontal stub_1_10 Fold_change" class="gt_row gt_right">-1.74</td>
<td headers="Frontal stub_1_10 Anova (p)" class="gt_row gt_right">0.02</td>
<td headers="Frontal stub_1_10 Reviewed" class="gt_row gt_left" style="font-size: 8;">unreviewed</td>
<td headers="Frontal stub_1_10 Protein.families" class="gt_row gt_left">V-ATPase H subunit family</td></tr>
    <tr><th id="stub_1_11" scope="row" class="gt_row gt_left gt_stub">Down</th>
<td headers="Frontal stub_1_11 Gene_name" class="gt_row gt_left" style="font-style: italic;">Bod1l</td>
<td headers="Frontal stub_1_11 Entry" class="gt_row gt_left">NA</td>
<td headers="Frontal stub_1_11 Uniprot_entry_name" class="gt_row gt_left">NA</td>
<td headers="Frontal stub_1_11 Full_name" class="gt_row gt_left">NA</td>
<td headers="Frontal stub_1_11 Accession" class="gt_row gt_left" style="font-size: 8;">BD1L1_MOUSE</td>
<td headers="Frontal stub_1_11 Description" class="gt_row gt_left" style="font-size: 8;">Biorientation of chromosomes in cell division protein 1-like 1 OS=Mus musculus GN=Bod1l PE=1 SV=1</td>
<td headers="Frontal stub_1_11 Fold_change" class="gt_row gt_right">-2.66</td>
<td headers="Frontal stub_1_11 Anova (p)" class="gt_row gt_right">0.03</td>
<td headers="Frontal stub_1_11 Reviewed" class="gt_row gt_left" style="font-size: 8;">NA</td>
<td headers="Frontal stub_1_11 Protein.families" class="gt_row gt_left">NA</td></tr>
    <tr><th id="stub_1_12" scope="row" class="gt_row gt_left gt_stub">Down</th>
<td headers="Frontal stub_1_12 Gene_name" class="gt_row gt_left" style="font-style: italic;">C1qbp</td>
<td headers="Frontal stub_1_12 Entry" class="gt_row gt_left">O35796</td>
<td headers="Frontal stub_1_12 Uniprot_entry_name" class="gt_row gt_left">C1QBP_RAT</td>
<td headers="Frontal stub_1_12 Full_name" class="gt_row gt_left">Complement component 1 Q subcomponent-binding protein, mitochondrial (GC1q-R protein) (Glycoprotein gC1qBP) (C1qBP)</td>
<td headers="Frontal stub_1_12 Accession" class="gt_row gt_left" style="font-size: 8;">C1QBP_MOUSE</td>
<td headers="Frontal stub_1_12 Description" class="gt_row gt_left" style="font-size: 8;">Complement component 1 Q subcomponent-binding protein, mitochondrial OS=Mus musculus GN=C1qbp PE=1 SV=1</td>
<td headers="Frontal stub_1_12 Fold_change" class="gt_row gt_right">-2.74</td>
<td headers="Frontal stub_1_12 Anova (p)" class="gt_row gt_right">0.05</td>
<td headers="Frontal stub_1_12 Reviewed" class="gt_row gt_left" style="font-size: 8;">reviewed</td>
<td headers="Frontal stub_1_12 Protein.families" class="gt_row gt_left">MAM33 family</td></tr>
    <tr><th id="stub_1_13" scope="row" class="gt_row gt_left gt_stub">Down</th>
<td headers="Frontal stub_1_13 Gene_name" class="gt_row gt_left" style="font-style: italic;">Capzb</td>
<td headers="Frontal stub_1_13 Entry" class="gt_row gt_left">Q5XI32</td>
<td headers="Frontal stub_1_13 Uniprot_entry_name" class="gt_row gt_left">CAPZB_RAT</td>
<td headers="Frontal stub_1_13 Full_name" class="gt_row gt_left">F-actin-capping protein subunit beta (CapZ beta)</td>
<td headers="Frontal stub_1_13 Accession" class="gt_row gt_left" style="font-size: 8;">CAPZB_RAT</td>
<td headers="Frontal stub_1_13 Description" class="gt_row gt_left" style="font-size: 8;">F-actin-capping protein subunit beta OS=Rattus norvegicus GN=Capzb PE=1 SV=1</td>
<td headers="Frontal stub_1_13 Fold_change" class="gt_row gt_right">-1.58</td>
<td headers="Frontal stub_1_13 Anova (p)" class="gt_row gt_right">0.04</td>
<td headers="Frontal stub_1_13 Reviewed" class="gt_row gt_left" style="font-size: 8;">reviewed</td>
<td headers="Frontal stub_1_13 Protein.families" class="gt_row gt_left">F-actin-capping protein beta subunit family</td></tr>
    <tr><th id="stub_1_14" scope="row" class="gt_row gt_left gt_stub">Down</th>
<td headers="Frontal stub_1_14 Gene_name" class="gt_row gt_left" style="font-style: italic;">Cct4</td>
<td headers="Frontal stub_1_14 Entry" class="gt_row gt_left">Q7TPB1</td>
<td headers="Frontal stub_1_14 Uniprot_entry_name" class="gt_row gt_left">TCPD_RAT</td>
<td headers="Frontal stub_1_14 Full_name" class="gt_row gt_left">T-complex protein 1 subunit delta (TCP-1-delta) (EC 3.6.1.-) (CCT-delta)</td>
<td headers="Frontal stub_1_14 Accession" class="gt_row gt_left" style="font-size: 8;">TCPD_RAT;TCPD_MOUSE</td>
<td headers="Frontal stub_1_14 Description" class="gt_row gt_left" style="font-size: 8;">T-complex protein 1 subunit delta OS=Rattus norvegicus GN=Cct4 PE=1 SV=3</td>
<td headers="Frontal stub_1_14 Fold_change" class="gt_row gt_right">-2.94</td>
<td headers="Frontal stub_1_14 Anova (p)" class="gt_row gt_right">0.05</td>
<td headers="Frontal stub_1_14 Reviewed" class="gt_row gt_left" style="font-size: 8;">reviewed</td>
<td headers="Frontal stub_1_14 Protein.families" class="gt_row gt_left">TCP-1 chaperonin family</td></tr>
    <tr><th id="stub_1_15" scope="row" class="gt_row gt_left gt_stub">Down</th>
<td headers="Frontal stub_1_15 Gene_name" class="gt_row gt_left" style="font-style: italic;">Cltb</td>
<td headers="Frontal stub_1_15 Entry" class="gt_row gt_left">P08082</td>
<td headers="Frontal stub_1_15 Uniprot_entry_name" class="gt_row gt_left">CLCB_RAT</td>
<td headers="Frontal stub_1_15 Full_name" class="gt_row gt_left">Clathrin light chain B (Lcb)</td>
<td headers="Frontal stub_1_15 Accession" class="gt_row gt_left" style="font-size: 8;">CLCB_MOUSE</td>
<td headers="Frontal stub_1_15 Description" class="gt_row gt_left" style="font-size: 8;">Clathrin light chain B OS=Mus musculus GN=Cltb PE=1 SV=1</td>
<td headers="Frontal stub_1_15 Fold_change" class="gt_row gt_right">-1.99</td>
<td headers="Frontal stub_1_15 Anova (p)" class="gt_row gt_right">0.05</td>
<td headers="Frontal stub_1_15 Reviewed" class="gt_row gt_left" style="font-size: 8;">reviewed</td>
<td headers="Frontal stub_1_15 Protein.families" class="gt_row gt_left">Clathrin light chain family</td></tr>
    <tr><th id="stub_1_16" scope="row" class="gt_row gt_left gt_stub">Down</th>
<td headers="Frontal stub_1_16 Gene_name" class="gt_row gt_left" style="font-style: italic;">Cs</td>
<td headers="Frontal stub_1_16 Entry" class="gt_row gt_left">Q8VHF5</td>
<td headers="Frontal stub_1_16 Uniprot_entry_name" class="gt_row gt_left">CISY_RAT</td>
<td headers="Frontal stub_1_16 Full_name" class="gt_row gt_left">Citrate synthase, mitochondrial (EC 2.3.3.1) (Citrate (Si)-synthase)</td>
<td headers="Frontal stub_1_16 Accession" class="gt_row gt_left" style="font-size: 8;">CISY_MOUSE</td>
<td headers="Frontal stub_1_16 Description" class="gt_row gt_left" style="font-size: 8;">Citrate synthase, mitochondrial OS=Mus musculus GN=Cs PE=1 SV=1</td>
<td headers="Frontal stub_1_16 Fold_change" class="gt_row gt_right">-2.01</td>
<td headers="Frontal stub_1_16 Anova (p)" class="gt_row gt_right">0.04</td>
<td headers="Frontal stub_1_16 Reviewed" class="gt_row gt_left" style="font-size: 8;">reviewed</td>
<td headers="Frontal stub_1_16 Protein.families" class="gt_row gt_left">Citrate synthase family</td></tr>
    <tr><th id="stub_1_17" scope="row" class="gt_row gt_left gt_stub">Down</th>
<td headers="Frontal stub_1_17 Gene_name" class="gt_row gt_left" style="font-style: italic;">Dnm3</td>
<td headers="Frontal stub_1_17 Entry" class="gt_row gt_left">Q08877</td>
<td headers="Frontal stub_1_17 Uniprot_entry_name" class="gt_row gt_left">DYN3_RAT</td>
<td headers="Frontal stub_1_17 Full_name" class="gt_row gt_left">Dynamin-3 (EC 3.6.5.5) (Dynamin, testicular) (T-dynamin)</td>
<td headers="Frontal stub_1_17 Accession" class="gt_row gt_left" style="font-size: 8;">DYN3_RAT;DYN3_MOUSE</td>
<td headers="Frontal stub_1_17 Description" class="gt_row gt_left" style="font-size: 8;">Dynamin-3 OS=Rattus norvegicus GN=Dnm3 PE=1 SV=2</td>
<td headers="Frontal stub_1_17 Fold_change" class="gt_row gt_right">-1.73</td>
<td headers="Frontal stub_1_17 Anova (p)" class="gt_row gt_right">0.01</td>
<td headers="Frontal stub_1_17 Reviewed" class="gt_row gt_left" style="font-size: 8;">reviewed</td>
<td headers="Frontal stub_1_17 Protein.families" class="gt_row gt_left">TRAFAC class dynamin-like GTPase superfamily, Dynamin/Fzo/YdjA family</td></tr>
    <tr><th id="stub_1_18" scope="row" class="gt_row gt_left gt_stub">Down</th>
<td headers="Frontal stub_1_18 Gene_name" class="gt_row gt_left" style="font-style: italic;">Gria2</td>
<td headers="Frontal stub_1_18 Entry" class="gt_row gt_left">P19491</td>
<td headers="Frontal stub_1_18 Uniprot_entry_name" class="gt_row gt_left">GRIA2_RAT</td>
<td headers="Frontal stub_1_18 Full_name" class="gt_row gt_left">Glutamate receptor 2 (GluR-2) (AMPA-selective glutamate receptor 2) (GluR-B) (GluR-K2) (Glutamate receptor ionotropic, AMPA 2)</td>
<td headers="Frontal stub_1_18 Accession" class="gt_row gt_left" style="font-size: 8;">GRIA2_RAT</td>
<td headers="Frontal stub_1_18 Description" class="gt_row gt_left" style="font-size: 8;">Glutamate receptor 2 OS=Rattus norvegicus GN=Gria2 PE=1 SV=2</td>
<td headers="Frontal stub_1_18 Fold_change" class="gt_row gt_right">-2.62</td>
<td headers="Frontal stub_1_18 Anova (p)" class="gt_row gt_right">0.03</td>
<td headers="Frontal stub_1_18 Reviewed" class="gt_row gt_left" style="font-size: 8;">reviewed</td>
<td headers="Frontal stub_1_18 Protein.families" class="gt_row gt_left">Glutamate-gated ion channel (TC 1.A.10.1) family, GRIA2 subfamily</td></tr>
    <tr><th id="stub_1_19" scope="row" class="gt_row gt_left gt_stub">Down</th>
<td headers="Frontal stub_1_19 Gene_name" class="gt_row gt_left" style="font-style: italic;">H2ax</td>
<td headers="Frontal stub_1_19 Entry" class="gt_row gt_left">A0ABK0LP06</td>
<td headers="Frontal stub_1_19 Uniprot_entry_name" class="gt_row gt_left">A0ABK0LP06_RAT</td>
<td headers="Frontal stub_1_19 Full_name" class="gt_row gt_left">H2A.X variant histone</td>
<td headers="Frontal stub_1_19 Accession" class="gt_row gt_left" style="font-size: 8;">H2AX_MOUSE</td>
<td headers="Frontal stub_1_19 Description" class="gt_row gt_left" style="font-size: 8;">Histone H2AX OS=Mus musculus GN=H2afx PE=1 SV=2</td>
<td headers="Frontal stub_1_19 Fold_change" class="gt_row gt_right">-2.46</td>
<td headers="Frontal stub_1_19 Anova (p)" class="gt_row gt_right">0.05</td>
<td headers="Frontal stub_1_19 Reviewed" class="gt_row gt_left" style="font-size: 8;">unreviewed</td>
<td headers="Frontal stub_1_19 Protein.families" class="gt_row gt_left">NA</td></tr>
    <tr><th id="stub_1_20" scope="row" class="gt_row gt_left gt_stub">Down</th>
<td headers="Frontal stub_1_20 Gene_name" class="gt_row gt_left" style="font-style: italic;">H3c1; H3c8; H3c10; H3c11</td>
<td headers="Frontal stub_1_20 Entry" class="gt_row gt_left">NA</td>
<td headers="Frontal stub_1_20 Uniprot_entry_name" class="gt_row gt_left">NA</td>
<td headers="Frontal stub_1_20 Full_name" class="gt_row gt_left">NA</td>
<td headers="Frontal stub_1_20 Accession" class="gt_row gt_left" style="font-size: 8;">H31_MOUSE;H3C_MOUSE</td>
<td headers="Frontal stub_1_20 Description" class="gt_row gt_left" style="font-size: 8;">Histone H3.1 OS=Mus musculus GN=Hist1h3a PE=1 SV=2</td>
<td headers="Frontal stub_1_20 Fold_change" class="gt_row gt_right">-2.74</td>
<td headers="Frontal stub_1_20 Anova (p)" class="gt_row gt_right">0.05</td>
<td headers="Frontal stub_1_20 Reviewed" class="gt_row gt_left" style="font-size: 8;">NA</td>
<td headers="Frontal stub_1_20 Protein.families" class="gt_row gt_left">NA</td></tr>
    <tr><th id="stub_1_21" scope="row" class="gt_row gt_left gt_stub">Down</th>
<td headers="Frontal stub_1_21 Gene_name" class="gt_row gt_left" style="font-style: italic;">Hadh</td>
<td headers="Frontal stub_1_21 Entry" class="gt_row gt_left">Q9WVK7</td>
<td headers="Frontal stub_1_21 Uniprot_entry_name" class="gt_row gt_left">HCDH_RAT</td>
<td headers="Frontal stub_1_21 Full_name" class="gt_row gt_left">Hydroxyacyl-coenzyme A dehydrogenase, mitochondrial (HCDH) (EC 1.1.1.35) (Medium and short-chain L-3-hydroxyacyl-coenzyme A dehydrogenase) (Short-chain 3-hydroxyacyl-CoA dehydrogenase)</td>
<td headers="Frontal stub_1_21 Accession" class="gt_row gt_left" style="font-size: 8;">HCDH_RAT</td>
<td headers="Frontal stub_1_21 Description" class="gt_row gt_left" style="font-size: 8;">Hydroxyacyl-coenzyme A dehydrogenase, mitochondrial OS=Rattus norvegicus GN=Hadh PE=2 SV=1</td>
<td headers="Frontal stub_1_21 Fold_change" class="gt_row gt_right">-2.63</td>
<td headers="Frontal stub_1_21 Anova (p)" class="gt_row gt_right">0.04</td>
<td headers="Frontal stub_1_21 Reviewed" class="gt_row gt_left" style="font-size: 8;">reviewed</td>
<td headers="Frontal stub_1_21 Protein.families" class="gt_row gt_left">3-hydroxyacyl-CoA dehydrogenase family</td></tr>
    <tr><th id="stub_1_22" scope="row" class="gt_row gt_left gt_stub">Down</th>
<td headers="Frontal stub_1_22 Gene_name" class="gt_row gt_left" style="font-style: italic;">Hadh</td>
<td headers="Frontal stub_1_22 Entry" class="gt_row gt_left">Q9WVK7</td>
<td headers="Frontal stub_1_22 Uniprot_entry_name" class="gt_row gt_left">HCDH_RAT</td>
<td headers="Frontal stub_1_22 Full_name" class="gt_row gt_left">Hydroxyacyl-coenzyme A dehydrogenase, mitochondrial (HCDH) (EC 1.1.1.35) (Medium and short-chain L-3-hydroxyacyl-coenzyme A dehydrogenase) (Short-chain 3-hydroxyacyl-CoA dehydrogenase)</td>
<td headers="Frontal stub_1_22 Accession" class="gt_row gt_left" style="font-size: 8;">HCDH_MOUSE</td>
<td headers="Frontal stub_1_22 Description" class="gt_row gt_left" style="font-size: 8;">Hydroxyacyl-coenzyme A dehydrogenase, mitochondrial OS=Mus musculus GN=Hadh PE=1 SV=2</td>
<td headers="Frontal stub_1_22 Fold_change" class="gt_row gt_right">-2.5</td>
<td headers="Frontal stub_1_22 Anova (p)" class="gt_row gt_right">0.05</td>
<td headers="Frontal stub_1_22 Reviewed" class="gt_row gt_left" style="font-size: 8;">reviewed</td>
<td headers="Frontal stub_1_22 Protein.families" class="gt_row gt_left">3-hydroxyacyl-CoA dehydrogenase family</td></tr>
    <tr><th id="stub_1_23" scope="row" class="gt_row gt_left gt_stub">Down</th>
<td headers="Frontal stub_1_23 Gene_name" class="gt_row gt_left" style="font-style: italic;">Hnrnpa1</td>
<td headers="Frontal stub_1_23 Entry" class="gt_row gt_left">P04256</td>
<td headers="Frontal stub_1_23 Uniprot_entry_name" class="gt_row gt_left">ROA1_RAT</td>
<td headers="Frontal stub_1_23 Full_name" class="gt_row gt_left">Heterogeneous nuclear ribonucleoprotein A1 (hnRNP A1) (Helix-destabilizing protein) (HDP) (Single-strand RNA-binding protein) (hnRNP core protein A1) [Cleaved into: Heterogeneous nuclear ribonucleoprotein A1, N-terminally processed]</td>
<td headers="Frontal stub_1_23 Accession" class="gt_row gt_left" style="font-size: 8;">ROA1_MOUSE</td>
<td headers="Frontal stub_1_23 Description" class="gt_row gt_left" style="font-size: 8;">Heterogeneous nuclear ribonucleoprotein A1 OS=Mus musculus GN=Hnrnpa1 PE=1 SV=2</td>
<td headers="Frontal stub_1_23 Fold_change" class="gt_row gt_right">-2.04</td>
<td headers="Frontal stub_1_23 Anova (p)" class="gt_row gt_right">0.03</td>
<td headers="Frontal stub_1_23 Reviewed" class="gt_row gt_left" style="font-size: 8;">reviewed</td>
<td headers="Frontal stub_1_23 Protein.families" class="gt_row gt_left">NA</td></tr>
    <tr><th id="stub_1_24" scope="row" class="gt_row gt_left gt_stub">Down</th>
<td headers="Frontal stub_1_24 Gene_name" class="gt_row gt_left" style="font-style: italic;">Hnrnpd</td>
<td headers="Frontal stub_1_24 Entry" class="gt_row gt_left">Q9JJ54</td>
<td headers="Frontal stub_1_24 Uniprot_entry_name" class="gt_row gt_left">HNRPD_RAT</td>
<td headers="Frontal stub_1_24 Full_name" class="gt_row gt_left">Heterogeneous nuclear ribonucleoprotein D0 (hnRNP D0) (AU-rich element RNA-binding protein 1)</td>
<td headers="Frontal stub_1_24 Accession" class="gt_row gt_left" style="font-size: 8;">HNRPD_RAT</td>
<td headers="Frontal stub_1_24 Description" class="gt_row gt_left" style="font-size: 8;">Heterogeneous nuclear ribonucleoprotein D0 OS=Rattus norvegicus GN=Hnrnpd PE=1 SV=2</td>
<td headers="Frontal stub_1_24 Fold_change" class="gt_row gt_right">-2.41</td>
<td headers="Frontal stub_1_24 Anova (p)" class="gt_row gt_right">0.03</td>
<td headers="Frontal stub_1_24 Reviewed" class="gt_row gt_left" style="font-size: 8;">reviewed</td>
<td headers="Frontal stub_1_24 Protein.families" class="gt_row gt_left">NA</td></tr>
    <tr><th id="stub_1_25" scope="row" class="gt_row gt_left gt_stub">Down</th>
<td headers="Frontal stub_1_25 Gene_name" class="gt_row gt_left" style="font-style: italic;">Ilf3</td>
<td headers="Frontal stub_1_25 Entry" class="gt_row gt_left">Q9JIL3</td>
<td headers="Frontal stub_1_25 Uniprot_entry_name" class="gt_row gt_left">ILF3_RAT</td>
<td headers="Frontal stub_1_25 Full_name" class="gt_row gt_left">Interleukin enhancer-binding factor 3</td>
<td headers="Frontal stub_1_25 Accession" class="gt_row gt_left" style="font-size: 8;">ILF3_MOUSE</td>
<td headers="Frontal stub_1_25 Description" class="gt_row gt_left" style="font-size: 8;">Interleukin enhancer-binding factor 3 OS=Mus musculus GN=Ilf3 PE=1 SV=2</td>
<td headers="Frontal stub_1_25 Fold_change" class="gt_row gt_right">-2.29</td>
<td headers="Frontal stub_1_25 Anova (p)" class="gt_row gt_right">0.04</td>
<td headers="Frontal stub_1_25 Reviewed" class="gt_row gt_left" style="font-size: 8;">reviewed</td>
<td headers="Frontal stub_1_25 Protein.families" class="gt_row gt_left">NA</td></tr>
    <tr><th id="stub_1_26" scope="row" class="gt_row gt_left gt_stub">Down</th>
<td headers="Frontal stub_1_26 Gene_name" class="gt_row gt_left" style="font-style: italic;">Inpp4a</td>
<td headers="Frontal stub_1_26 Entry" class="gt_row gt_left">Q62784</td>
<td headers="Frontal stub_1_26 Uniprot_entry_name" class="gt_row gt_left">INP4A_RAT</td>
<td headers="Frontal stub_1_26 Full_name" class="gt_row gt_left">Inositol polyphosphate-4-phosphatase type I A (Inositol polyphosphate 4-phosphatase type I) (Type I inositol 3,4-bisphosphate 4-phosphatase) (EC 3.1.3.66)</td>
<td headers="Frontal stub_1_26 Accession" class="gt_row gt_left" style="font-size: 8;">INP4A_RAT;INP4A_MOUSE</td>
<td headers="Frontal stub_1_26 Description" class="gt_row gt_left" style="font-size: 8;">Type I inositol 3,4-bisphosphate 4-phosphatase OS=Rattus norvegicus GN=Inpp4a PE=1 SV=1</td>
<td headers="Frontal stub_1_26 Fold_change" class="gt_row gt_right">-1.98</td>
<td headers="Frontal stub_1_26 Anova (p)" class="gt_row gt_right">0.04</td>
<td headers="Frontal stub_1_26 Reviewed" class="gt_row gt_left" style="font-size: 8;">reviewed</td>
<td headers="Frontal stub_1_26 Protein.families" class="gt_row gt_left">Inositol 3,4-bisphosphate 4-phosphatase family</td></tr>
    <tr><th id="stub_1_27" scope="row" class="gt_row gt_left gt_stub">Down</th>
<td headers="Frontal stub_1_27 Gene_name" class="gt_row gt_left" style="font-style: italic;">Mog</td>
<td headers="Frontal stub_1_27 Entry" class="gt_row gt_left">Q63345</td>
<td headers="Frontal stub_1_27 Uniprot_entry_name" class="gt_row gt_left">MOG_RAT</td>
<td headers="Frontal stub_1_27 Full_name" class="gt_row gt_left">Myelin-oligodendrocyte glycoprotein</td>
<td headers="Frontal stub_1_27 Accession" class="gt_row gt_left" style="font-size: 8;">MOG_MOUSE</td>
<td headers="Frontal stub_1_27 Description" class="gt_row gt_left" style="font-size: 8;">Myelin-oligodendrocyte glycoprotein OS=Mus musculus GN=Mog PE=1 SV=1</td>
<td headers="Frontal stub_1_27 Fold_change" class="gt_row gt_right">-3.05</td>
<td headers="Frontal stub_1_27 Anova (p)" class="gt_row gt_right">0.04</td>
<td headers="Frontal stub_1_27 Reviewed" class="gt_row gt_left" style="font-size: 8;">reviewed</td>
<td headers="Frontal stub_1_27 Protein.families" class="gt_row gt_left">Immunoglobulin superfamily, BTN/MOG family</td></tr>
    <tr><th id="stub_1_28" scope="row" class="gt_row gt_left gt_stub">Down</th>
<td headers="Frontal stub_1_28 Gene_name" class="gt_row gt_left" style="font-style: italic;">Ndufs6</td>
<td headers="Frontal stub_1_28 Entry" class="gt_row gt_left">P52504</td>
<td headers="Frontal stub_1_28 Uniprot_entry_name" class="gt_row gt_left">NDUS6_RAT</td>
<td headers="Frontal stub_1_28 Full_name" class="gt_row gt_left">NADH dehydrogenase [ubiquinone] iron-sulfur protein 6, mitochondrial (Complex I-13kD-A) (CI-13kD-A) (NADH-ubiquinone oxidoreductase 13 kDa-A subunit)</td>
<td headers="Frontal stub_1_28 Accession" class="gt_row gt_left" style="font-size: 8;">NDUS6_MOUSE</td>
<td headers="Frontal stub_1_28 Description" class="gt_row gt_left" style="font-size: 8;">NADH dehydrogenase [ubiquinone] iron-sulfur protein 6, mitochondrial OS=Mus musculus GN=Ndufs6 PE=1 SV=2</td>
<td headers="Frontal stub_1_28 Fold_change" class="gt_row gt_right">-2.35</td>
<td headers="Frontal stub_1_28 Anova (p)" class="gt_row gt_right">0.02</td>
<td headers="Frontal stub_1_28 Reviewed" class="gt_row gt_left" style="font-size: 8;">reviewed</td>
<td headers="Frontal stub_1_28 Protein.families" class="gt_row gt_left">Complex I NDUFS6 subunit family</td></tr>
    <tr><th id="stub_1_29" scope="row" class="gt_row gt_left gt_stub">Down</th>
<td headers="Frontal stub_1_29 Gene_name" class="gt_row gt_left" style="font-style: italic;">Ndufv1</td>
<td headers="Frontal stub_1_29 Entry" class="gt_row gt_left">Q5XIH3</td>
<td headers="Frontal stub_1_29 Uniprot_entry_name" class="gt_row gt_left">Q5XIH3_RAT</td>
<td headers="Frontal stub_1_29 Full_name" class="gt_row gt_left">NADH dehydrogenase [ubiquinone] flavoprotein 1, mitochondrial (EC 7.1.1.2)</td>
<td headers="Frontal stub_1_29 Accession" class="gt_row gt_left" style="font-size: 8;">NDUV1_MOUSE</td>
<td headers="Frontal stub_1_29 Description" class="gt_row gt_left" style="font-size: 8;">NADH dehydrogenase [ubiquinone] flavoprotein 1, mitochondrial OS=Mus musculus GN=Ndufv1 PE=1 SV=1</td>
<td headers="Frontal stub_1_29 Fold_change" class="gt_row gt_right">-4.16</td>
<td headers="Frontal stub_1_29 Anova (p)" class="gt_row gt_right">0.04</td>
<td headers="Frontal stub_1_29 Reviewed" class="gt_row gt_left" style="font-size: 8;">unreviewed</td>
<td headers="Frontal stub_1_29 Protein.families" class="gt_row gt_left">Complex I 51 kDa subunit family</td></tr>
    <tr><th id="stub_1_30" scope="row" class="gt_row gt_left gt_stub">Down</th>
<td headers="Frontal stub_1_30 Gene_name" class="gt_row gt_left" style="font-style: italic;">Omp</td>
<td headers="Frontal stub_1_30 Entry" class="gt_row gt_left">P08523</td>
<td headers="Frontal stub_1_30 Uniprot_entry_name" class="gt_row gt_left">OMP_RAT</td>
<td headers="Frontal stub_1_30 Full_name" class="gt_row gt_left">Olfactory marker protein (Olfactory neuronal-specific protein)</td>
<td headers="Frontal stub_1_30 Accession" class="gt_row gt_left" style="font-size: 8;">OMP_RAT</td>
<td headers="Frontal stub_1_30 Description" class="gt_row gt_left" style="font-size: 8;">Olfactory marker protein OS=Rattus norvegicus GN=Omp PE=1 SV=2</td>
<td headers="Frontal stub_1_30 Fold_change" class="gt_row gt_right">-3.62</td>
<td headers="Frontal stub_1_30 Anova (p)" class="gt_row gt_right">0.03</td>
<td headers="Frontal stub_1_30 Reviewed" class="gt_row gt_left" style="font-size: 8;">reviewed</td>
<td headers="Frontal stub_1_30 Protein.families" class="gt_row gt_left">Olfactory marker protein family</td></tr>
    <tr><th id="stub_1_31" scope="row" class="gt_row gt_left gt_stub">Down</th>
<td headers="Frontal stub_1_31 Gene_name" class="gt_row gt_left" style="font-style: italic;">Pgrmc1</td>
<td headers="Frontal stub_1_31 Entry" class="gt_row gt_left">P70580</td>
<td headers="Frontal stub_1_31 Uniprot_entry_name" class="gt_row gt_left">PGRC1_RAT</td>
<td headers="Frontal stub_1_31 Full_name" class="gt_row gt_left">Membrane-associated progesterone receptor component 1 (mPR) (25-DX) (Acidic 25 kDa protein) (Ventral midline antigen) (VEMA)</td>
<td headers="Frontal stub_1_31 Accession" class="gt_row gt_left" style="font-size: 8;">PGRC1_MOUSE</td>
<td headers="Frontal stub_1_31 Description" class="gt_row gt_left" style="font-size: 8;">Membrane-associated progesterone receptor component 1 OS=Mus musculus GN=Pgrmc1 PE=1 SV=4</td>
<td headers="Frontal stub_1_31 Fold_change" class="gt_row gt_right">-2.76</td>
<td headers="Frontal stub_1_31 Anova (p)" class="gt_row gt_right">0.02</td>
<td headers="Frontal stub_1_31 Reviewed" class="gt_row gt_left" style="font-size: 8;">reviewed</td>
<td headers="Frontal stub_1_31 Protein.families" class="gt_row gt_left">Cytochrome b5 family, MAPR subfamily</td></tr>
    <tr><th id="stub_1_32" scope="row" class="gt_row gt_left gt_stub">Down</th>
<td headers="Frontal stub_1_32 Gene_name" class="gt_row gt_left" style="font-style: italic;">Plp1</td>
<td headers="Frontal stub_1_32 Entry" class="gt_row gt_left">P60203</td>
<td headers="Frontal stub_1_32 Uniprot_entry_name" class="gt_row gt_left">MYPR_RAT</td>
<td headers="Frontal stub_1_32 Full_name" class="gt_row gt_left">Myelin proteolipid protein (PLP) (Lipophilin)</td>
<td headers="Frontal stub_1_32 Accession" class="gt_row gt_left" style="font-size: 8;">MYPR_MOUSE</td>
<td headers="Frontal stub_1_32 Description" class="gt_row gt_left" style="font-size: 8;">Myelin proteolipid protein OS=Mus musculus GN=Plp1 PE=1 SV=2</td>
<td headers="Frontal stub_1_32 Fold_change" class="gt_row gt_right">-2.62</td>
<td headers="Frontal stub_1_32 Anova (p)" class="gt_row gt_right">0.04</td>
<td headers="Frontal stub_1_32 Reviewed" class="gt_row gt_left" style="font-size: 8;">reviewed</td>
<td headers="Frontal stub_1_32 Protein.families" class="gt_row gt_left">Myelin proteolipid protein family</td></tr>
    <tr><th id="stub_1_33" scope="row" class="gt_row gt_left gt_stub">Down</th>
<td headers="Frontal stub_1_33 Gene_name" class="gt_row gt_left" style="font-style: italic;">Prkca</td>
<td headers="Frontal stub_1_33 Entry" class="gt_row gt_left">P05696</td>
<td headers="Frontal stub_1_33 Uniprot_entry_name" class="gt_row gt_left">KPCA_RAT</td>
<td headers="Frontal stub_1_33 Full_name" class="gt_row gt_left">Protein kinase C alpha type (PKC-A) (PKC-alpha) (EC 2.7.11.13)</td>
<td headers="Frontal stub_1_33 Accession" class="gt_row gt_left" style="font-size: 8;">KPCA_MOUSE;KPCA_RAT</td>
<td headers="Frontal stub_1_33 Description" class="gt_row gt_left" style="font-size: 8;">Protein kinase C alpha type OS=Mus musculus GN=Prkca PE=1 SV=3</td>
<td headers="Frontal stub_1_33 Fold_change" class="gt_row gt_right">-3.78</td>
<td headers="Frontal stub_1_33 Anova (p)" class="gt_row gt_right">0.00</td>
<td headers="Frontal stub_1_33 Reviewed" class="gt_row gt_left" style="font-size: 8;">reviewed</td>
<td headers="Frontal stub_1_33 Protein.families" class="gt_row gt_left">Protein kinase superfamily, AGC Ser/Thr protein kinase family, PKC subfamily</td></tr>
    <tr><th id="stub_1_34" scope="row" class="gt_row gt_left gt_stub">Down</th>
<td headers="Frontal stub_1_34 Gene_name" class="gt_row gt_left" style="font-style: italic;">Prkce</td>
<td headers="Frontal stub_1_34 Entry" class="gt_row gt_left">P09216</td>
<td headers="Frontal stub_1_34 Uniprot_entry_name" class="gt_row gt_left">KPCE_RAT</td>
<td headers="Frontal stub_1_34 Full_name" class="gt_row gt_left">Protein kinase C epsilon type (EC 2.7.11.13) (nPKC-epsilon)</td>
<td headers="Frontal stub_1_34 Accession" class="gt_row gt_left" style="font-size: 8;">KPCE_MOUSE</td>
<td headers="Frontal stub_1_34 Description" class="gt_row gt_left" style="font-size: 8;">Protein kinase C epsilon type OS=Mus musculus GN=Prkce PE=1 SV=1</td>
<td headers="Frontal stub_1_34 Fold_change" class="gt_row gt_right">-2.01</td>
<td headers="Frontal stub_1_34 Anova (p)" class="gt_row gt_right">0.04</td>
<td headers="Frontal stub_1_34 Reviewed" class="gt_row gt_left" style="font-size: 8;">reviewed</td>
<td headers="Frontal stub_1_34 Protein.families" class="gt_row gt_left">Protein kinase superfamily, AGC Ser/Thr protein kinase family, PKC subfamily</td></tr>
    <tr><th id="stub_1_35" scope="row" class="gt_row gt_left gt_stub">Down</th>
<td headers="Frontal stub_1_35 Gene_name" class="gt_row gt_left" style="font-style: italic;">Sdha</td>
<td headers="Frontal stub_1_35 Entry" class="gt_row gt_left">Q920L2</td>
<td headers="Frontal stub_1_35 Uniprot_entry_name" class="gt_row gt_left">SDHA_RAT</td>
<td headers="Frontal stub_1_35 Full_name" class="gt_row gt_left">Succinate dehydrogenase [ubiquinone] flavoprotein subunit, mitochondrial (EC 1.3.5.1) (Flavoprotein subunit of complex II) (Fp) (Malate dehydrogenase [quinone] flavoprotein subunit) (EC 1.1.5.-)</td>
<td headers="Frontal stub_1_35 Accession" class="gt_row gt_left" style="font-size: 8;">SDHA_MOUSE</td>
<td headers="Frontal stub_1_35 Description" class="gt_row gt_left" style="font-size: 8;">Succinate dehydrogenase [ubiquinone] flavoprotein subunit, mitochondrial OS=Mus musculus GN=Sdha PE=1 SV=1</td>
<td headers="Frontal stub_1_35 Fold_change" class="gt_row gt_right">-2.19</td>
<td headers="Frontal stub_1_35 Anova (p)" class="gt_row gt_right">0.02</td>
<td headers="Frontal stub_1_35 Reviewed" class="gt_row gt_left" style="font-size: 8;">reviewed</td>
<td headers="Frontal stub_1_35 Protein.families" class="gt_row gt_left">FAD-dependent oxidoreductase 2 family, FRD/SDH subfamily</td></tr>
    <tr><th id="stub_1_36" scope="row" class="gt_row gt_left gt_stub">Down</th>
<td headers="Frontal stub_1_36 Gene_name" class="gt_row gt_left" style="font-style: italic;">Septin8</td>
<td headers="Frontal stub_1_36 Entry" class="gt_row gt_left">B0BNF1</td>
<td headers="Frontal stub_1_36 Uniprot_entry_name" class="gt_row gt_left">SEPT8_RAT</td>
<td headers="Frontal stub_1_36 Full_name" class="gt_row gt_left">Septin-8</td>
<td headers="Frontal stub_1_36 Accession" class="gt_row gt_left" style="font-size: 8;">SEP11_MOUSE;SEP10_RAT;SEPT8_RAT</td>
<td headers="Frontal stub_1_36 Description" class="gt_row gt_left" style="font-size: 8;">Septin-11 OS=Mus musculus GN=Sept11 PE=1 SV=4</td>
<td headers="Frontal stub_1_36 Fold_change" class="gt_row gt_right">-1.94</td>
<td headers="Frontal stub_1_36 Anova (p)" class="gt_row gt_right">0.04</td>
<td headers="Frontal stub_1_36 Reviewed" class="gt_row gt_left" style="font-size: 8;">reviewed</td>
<td headers="Frontal stub_1_36 Protein.families" class="gt_row gt_left">TRAFAC class TrmE-Era-EngA-EngB-Septin-like GTPase superfamily, Septin GTPase family</td></tr>
    <tr><th id="stub_1_37" scope="row" class="gt_row gt_left gt_stub">Down</th>
<td headers="Frontal stub_1_37 Gene_name" class="gt_row gt_left" style="font-style: italic;">Serpina3m</td>
<td headers="Frontal stub_1_37 Entry" class="gt_row gt_left">Q63556</td>
<td headers="Frontal stub_1_37 Uniprot_entry_name" class="gt_row gt_left">SPA3M_RAT</td>
<td headers="Frontal stub_1_37 Full_name" class="gt_row gt_left">Serine protease inhibitor A3M (Serpin A3M) (Serine protease inhibitor 2.4) (SPI-2.4)</td>
<td headers="Frontal stub_1_37 Accession" class="gt_row gt_left" style="font-size: 8;">SPA3M_MOUSE</td>
<td headers="Frontal stub_1_37 Description" class="gt_row gt_left" style="font-size: 8;">Serine protease inhibitor A3M OS=Mus musculus GN=Serpina3m PE=1 SV=2</td>
<td headers="Frontal stub_1_37 Fold_change" class="gt_row gt_right">-2.58</td>
<td headers="Frontal stub_1_37 Anova (p)" class="gt_row gt_right">0.04</td>
<td headers="Frontal stub_1_37 Reviewed" class="gt_row gt_left" style="font-size: 8;">reviewed</td>
<td headers="Frontal stub_1_37 Protein.families" class="gt_row gt_left">Serpin family</td></tr>
    <tr><th id="stub_1_38" scope="row" class="gt_row gt_left gt_stub">Down</th>
<td headers="Frontal stub_1_38 Gene_name" class="gt_row gt_left" style="font-style: italic;">Slc8a2</td>
<td headers="Frontal stub_1_38 Entry" class="gt_row gt_left">P48768</td>
<td headers="Frontal stub_1_38 Uniprot_entry_name" class="gt_row gt_left">NAC2_RAT</td>
<td headers="Frontal stub_1_38 Full_name" class="gt_row gt_left">Sodium/calcium exchanger 2 (Na(+)/Ca(2+)-exchange protein 2) (Solute carrier family 8 member 2)</td>
<td headers="Frontal stub_1_38 Accession" class="gt_row gt_left" style="font-size: 8;">NAC2_RAT</td>
<td headers="Frontal stub_1_38 Description" class="gt_row gt_left" style="font-size: 8;">Sodium/calcium exchanger 2 OS=Rattus norvegicus GN=Slc8a2 PE=1 SV=1</td>
<td headers="Frontal stub_1_38 Fold_change" class="gt_row gt_right">-1.49</td>
<td headers="Frontal stub_1_38 Anova (p)" class="gt_row gt_right">0.05</td>
<td headers="Frontal stub_1_38 Reviewed" class="gt_row gt_left" style="font-size: 8;">reviewed</td>
<td headers="Frontal stub_1_38 Protein.families" class="gt_row gt_left">Ca(2+):cation antiporter (CaCA) (TC 2.A.19) family, SLC8 subfamily</td></tr>
    <tr><th id="stub_1_39" scope="row" class="gt_row gt_left gt_stub">Down</th>
<td headers="Frontal stub_1_39 Gene_name" class="gt_row gt_left" style="font-style: italic;">Stag3</td>
<td headers="Frontal stub_1_39 Entry" class="gt_row gt_left">Q99M76</td>
<td headers="Frontal stub_1_39 Uniprot_entry_name" class="gt_row gt_left">STAG3_RAT</td>
<td headers="Frontal stub_1_39 Full_name" class="gt_row gt_left">Cohesin subunit SA-3 (SCC3 homolog 3) (Stromal antigen 3) (Stromalin-3)</td>
<td headers="Frontal stub_1_39 Accession" class="gt_row gt_left" style="font-size: 8;">STAG3_RAT</td>
<td headers="Frontal stub_1_39 Description" class="gt_row gt_left" style="font-size: 8;">Cohesin subunit SA-3 OS=Rattus norvegicus GN=Stag3 PE=2 SV=1</td>
<td headers="Frontal stub_1_39 Fold_change" class="gt_row gt_right">-1.95</td>
<td headers="Frontal stub_1_39 Anova (p)" class="gt_row gt_right">0.00</td>
<td headers="Frontal stub_1_39 Reviewed" class="gt_row gt_left" style="font-size: 8;">reviewed</td>
<td headers="Frontal stub_1_39 Protein.families" class="gt_row gt_left">SCC3 family</td></tr>
    <tr><th id="stub_1_40" scope="row" class="gt_row gt_left gt_stub">Down</th>
<td headers="Frontal stub_1_40 Gene_name" class="gt_row gt_left" style="font-style: italic;">Stag3</td>
<td headers="Frontal stub_1_40 Entry" class="gt_row gt_left">Q99M76</td>
<td headers="Frontal stub_1_40 Uniprot_entry_name" class="gt_row gt_left">STAG3_RAT</td>
<td headers="Frontal stub_1_40 Full_name" class="gt_row gt_left">Cohesin subunit SA-3 (SCC3 homolog 3) (Stromal antigen 3) (Stromalin-3)</td>
<td headers="Frontal stub_1_40 Accession" class="gt_row gt_left" style="font-size: 8;">STAG3_MOUSE</td>
<td headers="Frontal stub_1_40 Description" class="gt_row gt_left" style="font-size: 8;">Cohesin subunit SA-3 OS=Mus musculus GN=Stag3 PE=1 SV=2</td>
<td headers="Frontal stub_1_40 Fold_change" class="gt_row gt_right">-1.81</td>
<td headers="Frontal stub_1_40 Anova (p)" class="gt_row gt_right">0.03</td>
<td headers="Frontal stub_1_40 Reviewed" class="gt_row gt_left" style="font-size: 8;">reviewed</td>
<td headers="Frontal stub_1_40 Protein.families" class="gt_row gt_left">SCC3 family</td></tr>
    <tr><th id="stub_1_41" scope="row" class="gt_row gt_left gt_stub">Down</th>
<td headers="Frontal stub_1_41 Gene_name" class="gt_row gt_left" style="font-style: italic;">Stx2</td>
<td headers="Frontal stub_1_41 Entry" class="gt_row gt_left">P50279</td>
<td headers="Frontal stub_1_41 Uniprot_entry_name" class="gt_row gt_left">STX2_RAT</td>
<td headers="Frontal stub_1_41 Full_name" class="gt_row gt_left">Syntaxin-2 (Epimorphin)</td>
<td headers="Frontal stub_1_41 Accession" class="gt_row gt_left" style="font-size: 8;">STX2_RAT</td>
<td headers="Frontal stub_1_41 Description" class="gt_row gt_left" style="font-size: 8;">Syntaxin-2 OS=Rattus norvegicus GN=Stx2 PE=1 SV=2</td>
<td headers="Frontal stub_1_41 Fold_change" class="gt_row gt_right">-4.53</td>
<td headers="Frontal stub_1_41 Anova (p)" class="gt_row gt_right">0.03</td>
<td headers="Frontal stub_1_41 Reviewed" class="gt_row gt_left" style="font-size: 8;">reviewed</td>
<td headers="Frontal stub_1_41 Protein.families" class="gt_row gt_left">Syntaxin family</td></tr>
    <tr><th id="stub_1_42" scope="row" class="gt_row gt_left gt_stub">Down</th>
<td headers="Frontal stub_1_42 Gene_name" class="gt_row gt_left" style="font-style: italic;">Uqcrc2</td>
<td headers="Frontal stub_1_42 Entry" class="gt_row gt_left">P32551</td>
<td headers="Frontal stub_1_42 Uniprot_entry_name" class="gt_row gt_left">QCR2_RAT</td>
<td headers="Frontal stub_1_42 Full_name" class="gt_row gt_left">Cytochrome b-c1 complex subunit 2, mitochondrial (Complex III subunit 2) (Core protein II) (Ubiquinol-cytochrome-c reductase complex core protein 2)</td>
<td headers="Frontal stub_1_42 Accession" class="gt_row gt_left" style="font-size: 8;">QCR2_RAT</td>
<td headers="Frontal stub_1_42 Description" class="gt_row gt_left" style="font-size: 8;">Cytochrome b-c1 complex subunit 2, mitochondrial OS=Rattus norvegicus GN=Uqcrc2 PE=1 SV=2</td>
<td headers="Frontal stub_1_42 Fold_change" class="gt_row gt_right">-2.75</td>
<td headers="Frontal stub_1_42 Anova (p)" class="gt_row gt_right">0.04</td>
<td headers="Frontal stub_1_42 Reviewed" class="gt_row gt_left" style="font-size: 8;">reviewed</td>
<td headers="Frontal stub_1_42 Protein.families" class="gt_row gt_left">Peptidase M16 family, UQCRC2/QCR2 subfamily</td></tr>
    <tr><th id="stub_1_43" scope="row" class="gt_row gt_left gt_stub">Down</th>
<td headers="Frontal stub_1_43 Gene_name" class="gt_row gt_left" style="font-style: italic;">Zfp11</td>
<td headers="Frontal stub_1_43 Entry" class="gt_row gt_left">A0A0G2K5Y5</td>
<td headers="Frontal stub_1_43 Uniprot_entry_name" class="gt_row gt_left">A0A0G2K5Y5_RAT</td>
<td headers="Frontal stub_1_43 Full_name" class="gt_row gt_left">Zinc finger protein 11</td>
<td headers="Frontal stub_1_43 Accession" class="gt_row gt_left" style="font-size: 8;">ZFP11_MOUSE</td>
<td headers="Frontal stub_1_43 Description" class="gt_row gt_left" style="font-size: 8;">Zinc finger protein 11 OS=Mus musculus GN=Zfp11 PE=2 SV=2</td>
<td headers="Frontal stub_1_43 Fold_change" class="gt_row gt_right">-2.52</td>
<td headers="Frontal stub_1_43 Anova (p)" class="gt_row gt_right">0.00</td>
<td headers="Frontal stub_1_43 Reviewed" class="gt_row gt_left" style="font-size: 8;">unreviewed</td>
<td headers="Frontal stub_1_43 Protein.families" class="gt_row gt_left">Krueppel C2H2-type zinc-finger protein family</td></tr>
    <tr><th id="stub_1_44" scope="row" class="gt_row gt_left gt_stub">Down</th>
<td headers="Frontal stub_1_44 Gene_name" class="gt_row gt_left" style="font-style: italic;">NA</td>
<td headers="Frontal stub_1_44 Entry" class="gt_row gt_left">Q6LED0</td>
<td headers="Frontal stub_1_44 Uniprot_entry_name" class="gt_row gt_left">H31_RAT</td>
<td headers="Frontal stub_1_44 Full_name" class="gt_row gt_left">Histone H3.1</td>
<td headers="Frontal stub_1_44 Accession" class="gt_row gt_left" style="font-size: 8;">CMC1_MOUSE</td>
<td headers="Frontal stub_1_44 Description" class="gt_row gt_left" style="font-size: 8;">Calcium-binding mitochondrial carrier protein Aralar1 OS=Mus musculus GN=Slc25a12 PE=1 SV=1</td>
<td headers="Frontal stub_1_44 Fold_change" class="gt_row gt_right">-1.79</td>
<td headers="Frontal stub_1_44 Anova (p)" class="gt_row gt_right">0.04</td>
<td headers="Frontal stub_1_44 Reviewed" class="gt_row gt_left" style="font-size: 8;">reviewed</td>
<td headers="Frontal stub_1_44 Protein.families" class="gt_row gt_left">Histone H3 family</td></tr>
    <tr class="gt_group_heading_row">
      <th colspan="11" class="gt_group_heading" scope="colgroup" id="Hippocampus">Hippocampus</th>
    </tr>
    <tr class="gt_row_group_first"><th id="stub_1_45" scope="row" class="gt_row gt_left gt_stub">Down</th>
<td headers="Hippocampus stub_1_45 Gene_name" class="gt_row gt_left" style="font-style: italic;">Atp5f1d</td>
<td headers="Hippocampus stub_1_45 Entry" class="gt_row gt_left">P35434</td>
<td headers="Hippocampus stub_1_45 Uniprot_entry_name" class="gt_row gt_left">ATPD_RAT</td>
<td headers="Hippocampus stub_1_45 Full_name" class="gt_row gt_left">ATP synthase F(1) complex subunit delta, mitochondrial (ATP synthase F1 subunit delta) (F-ATPase delta subunit)</td>
<td headers="Hippocampus stub_1_45 Accession" class="gt_row gt_left" style="font-size: 8;">ATPD_MOUSE</td>
<td headers="Hippocampus stub_1_45 Description" class="gt_row gt_left" style="font-size: 8;">ATP synthase subunit delta, mitochondrial OS=Mus musculus GN=Atp5d PE=1 SV=1</td>
<td headers="Hippocampus stub_1_45 Fold_change" class="gt_row gt_right">-2.15</td>
<td headers="Hippocampus stub_1_45 Anova (p)" class="gt_row gt_right">0.05</td>
<td headers="Hippocampus stub_1_45 Reviewed" class="gt_row gt_left" style="font-size: 8;">reviewed</td>
<td headers="Hippocampus stub_1_45 Protein.families" class="gt_row gt_left">ATPase epsilon chain family</td></tr>
    <tr><th id="stub_1_46" scope="row" class="gt_row gt_left gt_stub">Down</th>
<td headers="Hippocampus stub_1_46 Gene_name" class="gt_row gt_left" style="font-style: italic;">Cpne3</td>
<td headers="Hippocampus stub_1_46 Entry" class="gt_row gt_left">D3ZLA3</td>
<td headers="Hippocampus stub_1_46 Uniprot_entry_name" class="gt_row gt_left">D3ZLA3_RAT</td>
<td headers="Hippocampus stub_1_46 Full_name" class="gt_row gt_left">Copine-3 (Copine III)</td>
<td headers="Hippocampus stub_1_46 Accession" class="gt_row gt_left" style="font-size: 8;">CPNE3_MOUSE</td>
<td headers="Hippocampus stub_1_46 Description" class="gt_row gt_left" style="font-size: 8;">Copine-3 OS=Mus musculus GN=Cpne3 PE=1 SV=2</td>
<td headers="Hippocampus stub_1_46 Fold_change" class="gt_row gt_right">-2.71</td>
<td headers="Hippocampus stub_1_46 Anova (p)" class="gt_row gt_right">0.03</td>
<td headers="Hippocampus stub_1_46 Reviewed" class="gt_row gt_left" style="font-size: 8;">unreviewed</td>
<td headers="Hippocampus stub_1_46 Protein.families" class="gt_row gt_left">Copine family</td></tr>
    <tr><th id="stub_1_47" scope="row" class="gt_row gt_left gt_stub">Down</th>
<td headers="Hippocampus stub_1_47 Gene_name" class="gt_row gt_left" style="font-style: italic;">Prkca</td>
<td headers="Hippocampus stub_1_47 Entry" class="gt_row gt_left">P05696</td>
<td headers="Hippocampus stub_1_47 Uniprot_entry_name" class="gt_row gt_left">KPCA_RAT</td>
<td headers="Hippocampus stub_1_47 Full_name" class="gt_row gt_left">Protein kinase C alpha type (PKC-A) (PKC-alpha) (EC 2.7.11.13)</td>
<td headers="Hippocampus stub_1_47 Accession" class="gt_row gt_left" style="font-size: 8;">KPCA_MOUSE;KPCA_RAT</td>
<td headers="Hippocampus stub_1_47 Description" class="gt_row gt_left" style="font-size: 8;">Protein kinase C alpha type OS=Mus musculus GN=Prkca PE=1 SV=3</td>
<td headers="Hippocampus stub_1_47 Fold_change" class="gt_row gt_right">-3.68</td>
<td headers="Hippocampus stub_1_47 Anova (p)" class="gt_row gt_right">0.00</td>
<td headers="Hippocampus stub_1_47 Reviewed" class="gt_row gt_left" style="font-size: 8;">reviewed</td>
<td headers="Hippocampus stub_1_47 Protein.families" class="gt_row gt_left">Protein kinase superfamily, AGC Ser/Thr protein kinase family, PKC subfamily</td></tr>
    <tr><th id="stub_1_48" scope="row" class="gt_row gt_left gt_stub">Down</th>
<td headers="Hippocampus stub_1_48 Gene_name" class="gt_row gt_left" style="font-style: italic;">Ube2n</td>
<td headers="Hippocampus stub_1_48 Entry" class="gt_row gt_left">Q9EQX9</td>
<td headers="Hippocampus stub_1_48 Uniprot_entry_name" class="gt_row gt_left">UBE2N_RAT</td>
<td headers="Hippocampus stub_1_48 Full_name" class="gt_row gt_left">Ubiquitin-conjugating enzyme E2 N (EC 2.3.2.23) (Bendless-like ubiquitin-conjugating enzyme) (E2 ubiquitin-conjugating enzyme N) (Ubiquitin carrier protein N) (Ubiquitin-protein ligase N)</td>
<td headers="Hippocampus stub_1_48 Accession" class="gt_row gt_left" style="font-size: 8;">UBE2N_MOUSE</td>
<td headers="Hippocampus stub_1_48 Description" class="gt_row gt_left" style="font-size: 8;">Ubiquitin-conjugating enzyme E2 N OS=Mus musculus GN=Ube2n PE=1 SV=1</td>
<td headers="Hippocampus stub_1_48 Fold_change" class="gt_row gt_right">-2.94</td>
<td headers="Hippocampus stub_1_48 Anova (p)" class="gt_row gt_right">0.02</td>
<td headers="Hippocampus stub_1_48 Reviewed" class="gt_row gt_left" style="font-size: 8;">reviewed</td>
<td headers="Hippocampus stub_1_48 Protein.families" class="gt_row gt_left">Ubiquitin-conjugating enzyme family</td></tr>
    <tr><th id="stub_1_49" scope="row" class="gt_row gt_left gt_stub">Down</th>
<td headers="Hippocampus stub_1_49 Gene_name" class="gt_row gt_left" style="font-style: italic;">NA</td>
<td headers="Hippocampus stub_1_49 Entry" class="gt_row gt_left">Q6LED0</td>
<td headers="Hippocampus stub_1_49 Uniprot_entry_name" class="gt_row gt_left">H31_RAT</td>
<td headers="Hippocampus stub_1_49 Full_name" class="gt_row gt_left">Histone H3.1</td>
<td headers="Hippocampus stub_1_49 Accession" class="gt_row gt_left" style="font-size: 8;">CMC1_MOUSE</td>
<td headers="Hippocampus stub_1_49 Description" class="gt_row gt_left" style="font-size: 8;">Calcium-binding mitochondrial carrier protein Aralar1 OS=Mus musculus GN=Slc25a12 PE=1 SV=1</td>
<td headers="Hippocampus stub_1_49 Fold_change" class="gt_row gt_right">-1.6</td>
<td headers="Hippocampus stub_1_49 Anova (p)" class="gt_row gt_right">0.03</td>
<td headers="Hippocampus stub_1_49 Reviewed" class="gt_row gt_left" style="font-size: 8;">reviewed</td>
<td headers="Hippocampus stub_1_49 Protein.families" class="gt_row gt_left">Histone H3 family</td></tr>
  </tbody>
  &#10;</table>
</div>

``` r
  # opt_stylize(style = 1, color = 'gray', add_row_striping = TRUE)
```

## Plots

Quick Key word analysis \| based off the one from arroyo

``` r
key_annos <- read_excel("input/keywords_all_2025_12_09_uniprot.xlsx")
```

``` r
keywords_df <- lit_review_proteomics_sig %>%
  mutate(hashtags = str_split(Keywords, ";")) %>% # splits into list by sep ;
  tidyr::unnest(hashtags)
```

``` r
keywords_counts <- keywords_df %>%
  group_by(Region, In_EE) %>%
  filter(hashtags != "Reference proteome") %>%  # ignore ref
  dplyr::count(hashtags, sort = TRUE) %>%
  left_join(key_annos, by = join_by(hashtags == Name)) %>%
  mutate(In_EE = as.factor(In_EE),
         Category = as.factor(Category)) #%>%
```

``` r
summary(keywords_counts)
```

    ##     Region           In_EE       hashtags               n         
    ##  Length:175         Down:175   Length:175         Min.   : 1.000  
    ##  Class :character              Class :character   1st Qu.: 1.000  
    ##  Mode  :character              Mode  :character   Median : 2.000  
    ##                                                   Mean   : 3.029  
    ##                                                   3rd Qu.: 3.000  
    ##                                                   Max.   :28.000  
    ##                                                                   
    ##   Keyword ID                      Category  Gene Ontologies   
    ##  Length:175         Biological process:48   Length:175        
    ##  Class :character   Cellular component:34   Class :character  
    ##  Mode  :character   Molecular function:26   Mode  :character  
    ##                     Ligand            :24                     
    ##                     PTM               :22                     
    ##                     Domain            :11                     
    ##                     (Other)           :10                     
    ##   Definition          Synonyms            Links             Children        
    ##  Length:175         Length:175         Length:175         Length:175        
    ##  Class :character   Class :character   Class :character   Class :character  
    ##  Mode  :character   Mode  :character   Mode  :character   Mode  :character  
    ##                                                                             
    ##                                                                             
    ##                                                                             
    ##                                                                             
    ##    Parents           Statistics       
    ##  Length:175         Length:175        
    ##  Class :character   Class :character  
    ##  Mode  :character   Mode  :character  
    ##                                       
    ##                                       
    ##                                       
    ## 

``` r
plot(keywords_counts$n)
```

![](Shaw_2020_mice_proteomics_POSvsNEG_data_files/figure-gfm/unnamed-chunk-37-1.png)<!-- -->

``` r
keyplot <- keywords_counts %>%
  # ungroup() %>%
  # filter(percent > 2) %>% # 11 keywords
  filter(n > 2) %>%
  ggplot(aes(x = hashtags,
             # x = reorder(Category, n),
             y = n,
             # fill = reorder(hashtags, n))
             alpha = hashtags, 
             fill = Category)) +
  geom_col(color = "grey10", alpha = 0.9) +
  # facet_wrap(~In_EE, scales = "free_y", space = "free_y") +
  facet_nested(~Region+Category, 
             scales = 'free_x',
             space = 'free',
             # switch = "both"
             ) +
  labs(x = "Uniprot Keyword", y = "Count",
       title = "Shaw 2020 - Keywords in Significant Proteins",
       subtitle = "3 or more") +
  scale_fill_brewer(palette = "Spectral") +
  # # scale_fill_distiller(palette = "Spectral") +
  # scale_fill_paletteer("MetBrewer::Signac") +
  # coord_flip() +
  # theme_classic() +
  theme_bw() +
  theme(
    legend.position = "bottom",
     axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, # 90, hj1 t/R, vj.5 cent
                               size = 10, face = 'italic'), 
    axis.text.y = element_text(hjust = 1, color = "grey20", size = 10), # rota 90 & hjust=1 (top/right)
    )

keyplot
```

![](Shaw_2020_mice_proteomics_POSvsNEG_data_files/figure-gfm/unnamed-chunk-38-1.png)<!-- -->

------------------------------------------------------------------------

``` r
sessionInfo()
```

    ## R version 4.5.2 (2025-10-31)
    ## Platform: x86_64-pc-linux-gnu
    ## Running under: Ubuntu 22.04.5 LTS
    ## 
    ## Matrix products: default
    ## BLAS:   /usr/lib/x86_64-linux-gnu/blas/libblas.so.3.10.0 
    ## LAPACK: /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.10.0  LAPACK version 3.10.0
    ## 
    ## locale:
    ##  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
    ##  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
    ##  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
    ##  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
    ##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
    ## [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
    ## 
    ## time zone: Pacific/Auckland
    ## tzcode source: system (glibc)
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ##  [1] rmarkdown_2.30  readxl_1.4.5    ggh4x_0.3.1     pheatmap_1.0.13
    ##  [5] gt_1.2.0        ggrepel_0.9.6   lubridate_1.9.4 forcats_1.0.1  
    ##  [9] stringr_1.5.2   dplyr_1.1.4     purrr_1.1.0     readr_2.1.5    
    ## [13] tidyr_1.3.1     tibble_3.3.0    ggplot2_4.0.0   tidyverse_2.0.0
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] sass_0.4.10        utf8_1.2.6         generics_0.1.4     xml2_1.4.0        
    ##  [5] stringi_1.8.7      hms_1.1.3          digest_0.6.37      magrittr_2.0.4    
    ##  [9] evaluate_1.0.5     grid_4.5.2         timechange_0.3.0   RColorBrewer_1.1-3
    ## [13] fastmap_1.2.0      cellranger_1.1.0   scales_1.4.0       cli_3.6.5         
    ## [17] crayon_1.5.3       rlang_1.1.6        bit64_4.6.0-1      withr_3.0.2       
    ## [21] yaml_2.3.10        parallel_4.5.2     tools_4.5.2        tzdb_0.5.0        
    ## [25] vctrs_0.6.5        R6_2.6.1           lifecycle_1.0.4    fs_1.6.6          
    ## [29] bit_4.6.0          vroom_1.6.6        pkgconfig_2.0.3    pillar_1.11.1     
    ## [33] gtable_0.3.6       glue_1.8.0         Rcpp_1.1.0         xfun_0.53         
    ## [37] tidyselect_1.2.1   rstudioapi_0.17.1  knitr_1.50         farver_2.1.2      
    ## [41] htmltools_0.5.8.1  labeling_0.4.3     compiler_4.5.2     S7_0.2.1
