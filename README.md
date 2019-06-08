
<!-- README.md is generated from README.Rmd. Please edit that file -->
SAGA <img src="https://talbotsr.com/saga_package/logo.png" align="right" height="139" />
========================================================================================

### Surrogate Assay for Genotoxicity Assessment

<!-- badges: start -->
[![Travis build status](https://travis-ci.org/mytalbot/saga_package.svg?branch=master)](https://travis-ci.org/r-lib/usethis) <!-- badges: end -->

In our paper **"Predicting genotoxicity of integrating viral vectors for stem cell gene therapy using gene expression-based machine learning"** we report an improved in vitro test to determine the risk of insertional mutagenesis of integrating vectors for gene therapy. SAGA builds on the well accepted cell culture protocol of the in vitro immortalization assay (IVIM) but screens for the deregulation of oncogenic gene expression signatures. We demonstrate a new bioinformatic approach to correctly classify the mutagenic potential of retroviral vectors used in previous and current clinical trials.

Installation
------------

You can download the development version from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("mytalbot/saga_package")
library(saga)
```

### SAGA core data

You'll also need the SAGA core data. They are too large for the main package. Download them here (also from GitHub):

``` r
devtools::install_github("mytalbot/sagadata")
library(sagadata)
```

Example
-------

This is a basic example with just the data set. See the Vignette for a full example...

``` r
library(sagadata)
head(SAGA_Data)
#>       PROBE_ID X4467 X4469 X4471 X4473 X4907 X4908 X4909 X4990 X4991 X4992
#> 1 A_51_P100034  6484  6726  7105  6212  5632  6698  7014  3977  4214  5358
#> 2 A_51_P100174  7418  3014  6927  7533  2293  6385 15345  2853  2664  2824
#> 3 A_51_P100208    95    89    87    80    48    47    53    74    75    70
#> 4 A_51_P100289  5735  6423  6875  6702  4411  7013  6687  2473  2753  3580
#> 5 A_51_P100298  4999  2428  4851  3511  5192  6806  7089  2491  2276  2712
#> 6 A_51_P100309    78    87    72    70    44    39    42    79    80    86
#>   X4993 X4994 X4995 X4996 X4997 X5094 X5095 X5096 X5097 X5098 X5099 X5265
#> 1  5902  5202  4669  5060  5031  4545  4300  4789  4543  3997  4585  3639
#> 2  3577  3989  4057  4055  3623  6829  4505  9441  5114  3458  7110 10594
#> 3    92    70    68    66    72    60    73    60    54    52    51    67
#> 4  3861  3675  3358  3685  3434  3579  4366  4254  3544  4001  4362  4341
#> 5  3201  3643  3457  3314  3054  3594  2408  3575  3819  1962  3898  3288
#> 6    80    67    86    80    72    58    64    60    55    48    49    52
#>   X5268 X5271 X5399 X5400 X5401 X5402 X5403 X5404 X5405 X5406 X5407 X5408
#> 1  4628  5987  3205  3937  4310  4096  4312  3195  4923  4736  4153  3376
#> 2 33955 16809  9909 11927 10067 16787 11617 12640 12625 12270 11950 12741
#> 3    66    47    50    70    46    49    41    49    49    40    42    46
#> 4  8039  8340  5254  6806  6359  7594  6236  5825  8342  9862  6407  6280
#> 5  6668  3440  3800  2644  3090  3455  3280  2246  3785  2541  3550  2611
#> 6    58    46    48    59    52    53    46    50    52    52    44    42
#>   X5409 X5410 X5468 X5469 X5470 X5471 X5472 X5473 X5474 X5475 X5476 X5477
#> 1  4998  4694  4592  4441  4608  4996  3865  4712  3950  4689  5311  4604
#> 2 12917 14457 18165 12115 18627 13135 27971  8624 10236 24074  9166 30462
#> 3    49    48    47    50    47    44    40    41    48    48    40    45
#> 4  6997  7625 11971  8302  9580 12602 10240  9242  8001 10074 11064 10570
#> 5  3614  3120  3316  4093  4011  3702  2918  3388  3150  3586  3666  2904
#> 6    56    46    47    48    48    45    46    41    45    48    43    47
#>   X5478 X5479 X5602 X5603 X5604 X5605 X5606 X5607 X5608 X5609 X5610 X5611
#> 1  5563  4358  5121  4918  4439  4910  4564  5559  5760  6475  4706  4209
#> 2 31014 35166  6107  6836  8412  9944 12858 25444 20890 27830 28625 25634
#> 3    45    47    53    60    59    55    53    51    50    56    50    61
#> 4 14272 12015  4624  4085  4346  5304  4266  4282  3615  5154  3268  4196
#> 5  4121  3482  2745  2426  2820  3013  2014  3960  3276  3764  3196  2761
#> 6    44    44    62    55    54    53    47    84    49    51    52    53
#>   X5612 X5613 X5614 X5615 X5616 X5617 X5118 X5119 X5120 X5741 X5742 X5743
#> 1  5780  4953  5645  4861  5169  5484  3577  2912  3854  2376  2518  1800
#> 2 16214 15595 11898 15932 13151 35326  2467  2115  3572  2055  2394   972
#> 3    58    59    59    53    47    60    36    40    39    42    47    43
#> 4  3752  3749  4439  3180  3702  3708  2247  2101  2561  3875  3503  9958
#> 5  3089  1810  2359  2645  1561  2724  2503  2018  2841  2973  2923  1134
#> 6    52    54    62    55    48    53    39    41    41    60    56    53
#>   X5744 X5745 X5746 X5747 X5748 X9910 X9911 X9912 X9913 X9914 X9915 X9916
#> 1  2379  1843  1034  1790  2119  3198  4196  3369  4556  3599  4447  3665
#> 2  8789  8222  7033 12171  6320   766 21199  3412 26030 15125 20348  9341
#> 3    51    45    43    42    52    32    35    34    41    32    42    37
#> 4  3582  2706  2788  3284  4250  3716  5181  4222  5172  3715  4526  5355
#> 5  2312  2725  1794  3102  1650  1821  2422  1868  2171  2318  2105  2074
#> 6    99    79    48    72    91    37    41    32    38    32    42    36
#>   X9917 X6011 X6012 X6013 X6014 X6015 X6016 X6017 X6018
#> 1  3927  3506  3918  4200  3708  3800  3601  4038  3721
#> 2  4464 19612 15225 27329 12677 13359 14193 23387 21222
#> 3    38    49    41    38    41    48    37    38    38
#> 4  4286 11850 10543 13496  7624  6764  5556 10405  7397
#> 5  2045  3420  3295  4028  3158  3094  2766  3547  2403
#> 6    32    59    43    45    57    42    44    43    37
```
