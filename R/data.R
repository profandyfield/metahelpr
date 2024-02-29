#' Brewin & Field (2024) without effect sizes
#'
#' Data from a meta-analysis looking at predictors of disorganization in trauma memories. See also [brewin_full].
#'
#' Background: Taylor et al. (Clinical Psychological Science, 2022) reported that in healthy participants memories of
#' traumatic and comparison films did not differ in coherence. The lack of a group diagnosed with posttraumatic stress disorder (PTSD),
#' as well as limitations of the trauma film paradigm, mean that their design is unable to directly test predictions made by clinical theories of PTSD. Brewin and Field argued that
#' there is convincing evidence for trauma memories in PTSD being incoherent or disorganized and conducted a meta-analysis to estimate the effect size for PTSD status and memory incoherence/disorganization.
#' These are the raw data from this meta-analysis. Key to the analysis is the fact that there are two distinct approaches to measuring disorganization/incoherence with
#' acceptable face validity. One, based on therapy with PTSD patients, is to elicit a very detailed trauma narrative that includes the worst moments and have judges rate individual utterance units for markers of disorganization such as repetition and non-consecutive chunks
#' (The FOA method). The other approach, employed by Taylor et al. (2022) among others, has judges or participants rate the entire memory or narrative. The tibble contains the following variables
#'
#'   * **author**: the author of the paper
#'   * **study**: numeric study identifier
#'   * **measure**: the measure of disorganization used
#'   * **high_disorganization**: is high disorganization represented by a low score or a high score
#'   * **rater**: were memories rated by the participant (self) or researcher (judge)
#'   * **method**: the method for classifying memories
#'   * **foa_type**: if the FOA method was employed was it global or detailed?
#'   * **foa**:	was the FOA method used to measure disorganization?
#'   * **clinical_group**: Did the clinical group have diagnosis of PTSD or ASD
#'   * **age**:	Adult or youth sample?
#'   * **clinical_mean**:	mean disorganization in the clinical group
#'   * **clinical_sd**:	standard deviation of disorganization in the clinical group
#'   * **clinical_n**:	sample size disorganization in the clinical group
#'   * **control_mean**: mean disorganization in the control group
#'   * **control_sd**:	standard deviation of disorganization in the control group
#'   * **control_n**:	sample size disorganization in the control group
#'   * **es_id**:	numeric identifier for unique effect sizes
#'   * **t**:	*t*-statistic for the comparison of clinical and control groups on disorganisation scores
#'   * **r**:	correlation coefficient for PTSD/ASD measures and disorganization
#'   * **n_for_r**:	sample size for the correlation coefficient
#'   * **ptsd_measure**:	The PTSD measure used
#'   * **ptsd_mean_clin**: mean PTSD score on the PTSD measure in the clinical group
#'   * **ptsd_sd_clin**: standard deviation of the PTSD score on the PTSD measure in the clinical group
#'   * **ptsd_mean_ctrl**: mean PTSD score on the PTSD measure in the control group
#'   * **ptsd_sd_ctrl**: standard deviation of the PTSD score on the PTSD measure in the control group
#'   * **ptsd_mean_all**:	mean PTSD score on the PTSD measure for the whole sample
#'   * **ptsd_sd_all**:	standard deviation of the PTSD score on the PTSD measure for the whole sample
#'   * **source**: text describing source of the descriptive statistics in the paper
#'   * **notes**:	text notes about the study
#'   * **direction**:	the direction of the effect size
#'   * **foa_gp**: a categorical variable combining information about the `foa_type` and `foa` group.
#'   * **global_vs_detail**: the first of three dummy variables. This one compares global FOA methodology to detailed FOA methodology
#'   * **dis_vs_detail**:	the second of three dummy variables. This one compares non-FOA methodology that measures disorganization to detailed FOA methodology
#'   * **org_vs_detail**: the last of three dummy variables. This one compares non-FOA methodology that measures organization to detailed FOA methodology
#'
#' @docType data
#' @format A tibble with 80 rows and 34 variables.
#' @source [https://osf.io/597hr/](https://osf.io/597hr/)
#' @references
#'
#'  * Brewin, C. & Field, A. P. (submitted). Meta-analysis shows trauma memories in PTSD lack coherence: A response to Taylor et al. (2022). [https://osf.io/597hr/](https://osf.io/597hr/)
#'  * Taylor, A., Zajac, R., Takarangi, M. K. T., & Garry, M. (2022). Evidence from the trauma-film paradigm that traumatic and nontraumatic memories are statistically equivalent on coherence. *Clinical Psychological Science*, 10(3), 417–429. [https://doi.org/10.1177/21677026211053312](https://doi.org/10.1177/21677026211053312)
#'
"brewin_2024"

#' Brewin & Field (2024) with effect sizes
#'
#' Data from a meta-analysis looking at predictors of disorganization in trauma memories. This is a version of the data in [brewin_2024]
#' where the effect sizes have already been computed and some non-essential variables are omitted.
#'
#' Background: Taylor et al. (Clinical Psychological Science, 2022) reported that in healthy participants memories of
#' traumatic and comparison films did not differ in coherence. The lack of a group diagnosed with posttraumatic stress disorder (PTSD),
#' as well as limitations of the trauma film paradigm, mean that their design is unable to directly test predictions made by clinical theories of PTSD. Brewin and Field argued that
#' there is convincing evidence for trauma memories in PTSD being incoherent or disorganized and conducted a meta-analysis to estimate the effect size for PTSD status and memory incoherence/disorganization.
#' These are the raw data from this meta-analysis. Key to the analysis is the fact that there are two distinct approaches to measuring disorganization/incoherence with
#' acceptable face validity. One, based on therapy with PTSD patients, is to elicit a very detailed trauma narrative that includes the worst moments and have judges rate individual utterance units for markers of disorganization such as repetition and non-consecutive chunks
#' (The FOA method). The other approach, employed by Taylor et al. (2022) among others, has judges or participants rate the entire memory or narrative. The tibble contains the following variables
#'
#'   * **author**: the author of the paper
#'   * **study**: numeric study identifier
#'   * **high_disorganization**: is high disorganization represented by a low score or a high score
#'   * **clinical_group**: Did the clinical group have diagnosis of PTSD or ASD
#'   * **age**:	Adult or youth sample?
#'   * **clinical_mean**:	mean disorganization in the clinical group
#'   * **clinical_sd**:	standard deviation of disorganization in the clinical group
#'   * **clinical_n**:	sample size disorganization in the clinical group
#'   * **control_mean**: mean disorganization in the control group
#'   * **control_sd**:	standard deviation of disorganization in the control group
#'   * **control_n**:	sample size disorganization in the control group
#'   * **es_id**:	numeric identifier for unique effect sizes
#'   * **foa_gp**: a categorical variable combining information about the `foa_type` and `foa` group.
#'   * **global_vs_detail**: the first of three dummy variables. This one compares global FOA methodology to detailed FOA methodology
#'   * **dis_vs_detail**:	the second of three dummy variables. This one compares non-FOA methodology that measures disorganization to detailed FOA methodology
#'   * **org_vs_detail**: the last of three dummy variables. This one compares non-FOA methodology that measures organization to detailed FOA methodology
#'   * **g**: Hedge's g
#'   * **v_g**: variance for Hedge's g
#'   * **n_total**: total sample size
#'
#' @docType data
#' @format A tibble with 80 rows and 19 variables.
#' @source [https://osf.io/597hr/](https://osf.io/597hr/)
#' @references
#'
#'  * Brewin, C. & Field, A. P. (submitted). Meta-analysis shows trauma memories in PTSD lack coherence: A response to Taylor et al. (2022). [https://osf.io/597hr/](https://osf.io/597hr/)
#'  * Taylor, A., Zajac, R., Takarangi, M. K. T., & Garry, M. (2022). Evidence from the trauma-film paradigm that traumatic and nontraumatic memories are statistically equivalent on coherence. *Clinical Psychological Science*, 10(3), 417–429. [https://doi.org/10.1177/21677026211053312](https://doi.org/10.1177/21677026211053312)
#'
"brewin_es"

#' Pearce & Field (2016)
#'
#' Background: Despite a general perception that violent or scary television creates anxiety in children, the research literature was small and disparate.
#' In 2016 we conducted a meta-analysis that quantified the impact of scary television and film on children's internalizing emotions (fear, anxiety, sadness, and sleep problems).
#' `pearce_2016` contains the key variables from the data used in this study. The data contains the following variables
#'
#'   * **author**: the author of the paper
#'   * **study_id**: numeric study identifier
#'   * **es_id**:	numeric identifier for unique effect sizes
#'   * **year**: year of publication
#'   * **methodology**: Did the clinical group have diagnosis of PTSD or ASD
#'   * **responder**: was the responder the child, the parent or both
#'   * **age_range_cat**: age range of participants (text)
#'   * **age_range_num**:	age range of participants (numeric)
#'   * **mean_age_est**:	estimated mean age of participants
#'   * **all_under_10**: was all of the sample under 10 years old?
#'   * **measure_type**: Type of measure used (e.g. self-report)
#'   * **outcome**: Outcome measured (fear, PTSD, general, Depression/sadness, physiology, sleep disorder)
#'   * **news_fiction**: was the media consumed news or fiction?
#'   * **r**: effect size *r*
#'   * **var_r**:	variance of the effect size
#'   * **n**:	total sample size onw hich the effect size is based
#'   * **tv**: was the media consumed television or mixed media?
#'   * **content**: classification of the content consumed (e.g. terrorism)
#'   * **violence**: rating of how violent the content was from 0 (not at all violent) to 3 (most violent)
#'
#'
#' @docType data
#' @format A tibble with 129 rows and 19 variables.
#' @source [https://doi.org/doi:10.1111/hcre.12069](https://doi.org/doi:10.1111/hcre.12069)
#' @references
#'
#'  * Pearce, L. J., & Field, A. P. (2016). The impact of 'scary' TV and film on children’s internalizing emotions: A meta-analysis. *Human Communication Research*, 42, 98–121. [https://doi.org/doi:10.1111/hcre.12069](https://doi.org/doi:10.1111/hcre.12069)
#'
"pearce_2016"

