
collate_models <- function(model, digits = 3, p_digits = 3, row = 1, mv = T){
  fmt <- paste0("%.", digits, "f")
  
  out <- tidy(model, conf.int = TRUE)[row, ]
  
  if(mv){
    out <- out |>
      mutate(
        sigma_b = model$sigma2[1],
        sigma_w = model$sigma2[2])
  } else {
    out <- out |>
      mutate(
        tau_2 = model$tau2)
  }
  
  out <- out |>
    mutate(
      q = model$QE,
      q_p = metafor::fmtp(model$QEp, p_digits),
      across(where(is.numeric), \(x) sprintf(fmt = fmt, x)),
      ci = paste0("[", conf.low, ", ", conf.high, "]"),
      k = model$k,
      q_df = ifelse(mv, model$QEdf, NA),
    )
  
  if(mv){
    out |> 
      select(k, sigma_b, sigma_w, q, q_df, q_p, estimate, ci, statistic, p.value)
  } else {
    out |> 
      select(k, tau_2, q, q_p, estimate, ci, statistic, p.value)
  }
}

#' Fit individual meta-analyses for each level of a categorical predictor.
#' 
#' `get_mas()` fits individual meta-analyses for each level of a categorical variable and (optionally) collates the results into a tabulated form for printing in a quarto document.
#' Models are fitted using using the `rma.mv()` function from the `metafor` package.
#'
#' @param tibble A tibble containing the effect sizes, their variances,a study ID and an effect size ID
#' @param predictor Name of the categorical variable for which you want individual meta-analyses
#' @param es_id  Name of the variable that identifies unique effect sizes
#' @param es_id  Name of the variable that identifies unique studies
#' @param digits  number of decimal places to print in the output
#' @param p_digits  number of decimal places for *p*-values
#' @param summary if TRUE output a summary table, otherwise store the raw data and models
#' 
#' @return If `summary = TRUE` the function returns a tibble containing the results of the meta-analysis for each category in a separate row. 
#' Columns show 
#' 
#' * *k*: number of effect sizes
#' * \eqn{\sigma_b}: between study variance
#' * \eqn{\sigma_w}: within-study variance
#' * *q*: Q statistic for heterogeneity,
#' * \eqn{q_{df}}: degrees of freedom for the *Q* statistic,
#' * \eqn{q_p}: *p*-value for the Q statistic
#' * estimate: The pooled effect size estimate
#' * ci: the 95% confidence interval for the pooled effect size estimate
#' * statistic: the test statistic for the pooled effect size estimate
#' * p.value: the test statistic for the pooled effect size estimate
#'
#' If `summary = FALSE` the function returns a tibble containing the models. Each row represents a category of the predictor, and columns contain
#' 
#' * `name_of_predictor`: This column will be named after the predictor entered into the function and will list the categories (factor levels) in each row.
#' * data: contains the raw data on which the model is based (that is, the data for the particular category)
#' * model: each row contains the `rma.mv` object for the model fitted within each category
#' * coefs: each row contains a table of coefficients for the model fitted within each category
#' 

get_mas <- function(tibble, predictor, es_id, study_id, digits = 2, p_digits = 3, summary = TRUE){
  mas <- tibble  |>
    dplyr::arrange({{predictor}}) |> 
    dplyr::group_by({{predictor}})  |>  
    tidyr::nest()  |> 
    dplyr::mutate(
      model = purrr::map(.x = data, 
                         .f = \(es_tib) rma.mv(yi = g, V = v_g, random = ~1|{{study_id}}/{{es_id}}, data = es_tib)),
      coefs = purrr::map(.x = model, .f = \(m) collate_models(m, digits = digits, p_digits = p_digits))
    )
  
  if(summary){
    mas  |>
      dplyr::select(c({{predictor}}, coefs)) |> 
      tidyr::unnest(coefs)
  } else {
    mas
  }
}




#' Run individual meta-analyses with a moderator variable for each level of a categorical predictor.
#' 
#' `get_mod_mas()` does the same thing as [`get_mas()`] but allows the user to specify a moderator/predictor within the individual meta-analyses.
#' So it fits individual meta-analyses models with a moderator specified using the `rma.mv()` function from the `metafor` package 
#' within each level of a categorical variable. Optionally, the results can be collated into a tabulated form for printing in a quarto document. 
#' 
#' @param tibble A tibble containing the effect sizes, their variances,a study ID and an effect size ID
#' @param predictor Name of the categorical variable for which you want individual meta-analyses
#' @param moderator Name of the moderator variable to be included within each sub-analysis
#' @param es_id  Name of the variable that identifies unique effect sizes
#' @param es_id  Name of the variable that identifies unique studies
#' @param digits  number of decimal places to print in the output
#' @param p_digits  number of decimal places for *p*-values
#' @param summary if TRUE output a summary table, otherwise store the raw data and models
#' 
#' @return A tibble containing the results of the meta-analysis for each category in a separate row. 
#' Columns show 
#' 
#' * *k*: number of effect sizes
#' * \eqn{\sigma_b}: between study variance
#' * \eqn{\sigma_w}: within-study variance
#' * *q*: Q statistic for heterogeneity,
#' * \eqn{q_{df}}: degrees of freedom for the *Q* statistic,
#' * \eqn{q_p}: *p*-value for the Q statistic
#' * estimate: The estimate of the effect of the moderator
#' * ci: the 95% confidence interval for the estimate of the effect of the moderator
#' * statistic: the test statistic for the estimate of the effect of the moderator
#' * p.value: the test statistic for the estimate of the effect of the moderator
#'
#'

get_mod_mas <- function(tibble, predictor, moderator, es_id, study_id, digits = 2, p_digits = 3, summary = TRUE){
  mas <- tibble  |>
    dplyr::arrange({{predictor}}) |> 
    dplyr::group_by({{predictor}})  |>  
    tidyr::nest()  |> 
    dplyr::mutate(
      model = purrr::map(.x = data, 
                         .f = \(es_tib) rma.mv(yi = g, V = v_g, mods = ~{{moderator}}, random = ~1|{{study_id}}/{{es_id}}, data = es_tib)),
      coefs = purrr::map(.x = model, .f = \(m) collate_models(m, row = 2, digits = digits, p_digits = p_digits))
    )
  
  if(summary){
    mas  |>
      dplyr::select(c({{predictor}}, coefs)) |> 
      tidyr::unnest(coefs)
  } else {
    mas
  }
}




#' Report residual heterogeneity statistics
#' 
#' `report_het()` is a helper function that collates information from heterogeneity tests and outputs text that summarizes the results in
#' a format that will render in quarto. Basically it helps you to report results.
#' 
#' @param metaobject A meta-analyses object (`rma` or `rma.mv`) created using the `metafor` package
#' @param digits  number of decimal places to print in the output
#' @param p_digits  number of decimal places for *p*-values
#' 
#' @return A text string including LaTeX formatting for mathematical symbols that reports the *Q* statistic, its degrees of freedom and it's associated *p*-value.
#' 
#'
#'

report_het <- function(metaobject, digits = 2, p_digits = 3){
  paste0("$Q_E$(", metaobject$QEdf, ") = ", metafor::fmtx(metaobject$QE, digits), ", *p* ", metafor::fmtp(metaobject$QEp, p_digits, sep = T, equal = T))
}


# get row of tidy output from metafor

get_row <- function(tidyobject, row = 2, digits = 2, p_digits = 3){
  n <- nrow(tidyobject)
  
  tidyobject |>
    dplyr::mutate(
      row_number = 1:n,
      p.value = metafor::fmtp(p.value, p_digits),
      ci = paste0("[", metafor::fmtx(conf.low, digits), ", ", metafor::fmtx(conf.high, digits), "]"),
      across(.cols = where(is.double), \(x) metafor::fmtx(x, digits))
    ) |>
    dplyr::filter(row_number == row)
}


#' Report residual heterogeneity statistics
#' 
#' `report_mod()` is a helper function that outputs text that reports the omnibus statistical tests from a moderation model.
#' Basically it helps you to report results in quarto.
#' 
#' @param metaobject A meta-analyses object (`rma` or `rma.mv`) created using the `metafor` package
#' @param digits  number of decimal places to print in the output
#' @param p_digits  number of decimal places for *p*-values
#' 
#' @return A text string including LaTeX formatting for mathematical symbols that reports the *Q* statistic (or *F* for robust models), the associated degrees of freedom and *p*-value.
#' 


report_mod <- function(metaobject, digits = 2, p_digits = 3){
  
  if(is.null(metaobject$robumethod)){
    paste0("$Q_E$(", metaobject$QMdf[1], ") = ", metafor::fmtx(metaobject$QM, digits), ", *p* ", metafor::fmtp(metaobject$QMp, p_digits, sep = T, equal = T))
  } else {
    paste0("$F$(", metaobject$QMdf[1], ", ", metaobject$QMdf[2], ") = ", metafor::fmtx(metaobject$QM, digits), ", *p* ", metafor::fmtp(metaobject$QMp, p_digits, sep = T, equal = T))
  }
}

#' Report parameters from a meta-analysis object
#' 
#' `report_pars()` is a helper function that outputs text that reports the individual effects from a meta-analysis model.
#' Basically it helps you to report results in quarto.
#' 
#' @param metaobject A meta-analyses object (`rma` or `rma.mv`) created using the `metafor` package
#' @param row The row of output you wish to report (if you have more than one effect)
#' @param digits  number of decimal places to print in the output
#' @param p_digits  number of decimal places for *p*-values
#' @param mod Are you reporting the effect of a moderator variable? This affects the symbols used for the parameter estimates. For moderators the 
#' symbol \eqn{\hat{\beta}} is used, for overall effects \eqn{\hat{\theta}} is used.
#' 
#' @return A text string including LaTeX formatting for mathematical symbols that reports the estimate (\eqn{\hat{\beta}} for moderators, \eqn{\hat{\theta}} for overall effects),
#' 95% confidence interval, associated test-statistic and *p*-value.
#' 


report_pars <- function(metaobject, row = 1, digits = 2, p_digits = 3, mod = F){
  tidyobject <- tidy(metaobject, conf.int = T)
  row <- get_row(tidyobject, row, digits, p_digits)
  
  symbol <- ifelse(mod == TRUE, 
                 paste0("$\\hat{\\beta}$ = "),
                 paste0("$\\hat{\\theta}$ = "))
  
  stat <- ifelse(!is.null(metaobject$robumethod), 
                 paste0(", *t*(", metaobject$dfs, ") = "),
                 paste0(", *z* = "))

  
  paste0(symbol, row$estimate, " ", row$ci, stat, row$statistic, ", *p* ", metafor::fmtp(row$p.value, p_digits, sep = T, equal = T))
}

#' Tabulate parameters from a meta-analysis object
#' 
#' `report_par_tbl()` is a helper function that outputs a tibble of the table of coefficients of a meta-analyses object (`rma` or `rma.mv`) 
#' created using the `metafor` package. This will mostly be useful for models containing predictors of effect sizes (so-called meta-regression).
#' 
#' @param metaobject A meta-analyses object (`rma` or `rma.mv`) created using the `metafor` package
#' @param digits  number of decimal places to print in the output
#' @param p_digits  number of decimal places for *p*-values
#' @param pred_labels A vector of text labels to use for each row of the table. For example, `pred_labels = c("Intercept", "FOA (Global) vs. FOA (Detailed)", "Not FOA (Disorganisation) vs. FOA (Detailed)", "Not FOA (Organisation) vs. FOA (Detailed)")`
#' 
#' @return A tibble containing (for each predictor in the model) the name of the predictor (obtained from `pred_labels`), the effect size estimate, 95% confidence interval, test statistic and *p*-value.
#'


report_par_tbl <- function(metaobject, digits = 2, p_digits = 3, pred_labels){
  fmt <- paste0("%.", digits, "f")
  
  broom::tidy(metaobject, conf.int = T) |>
    mutate(
      p = metafor::fmtp(p.value, 3),
      across(where(is.numeric), \(x) sprintf(fmt = fmt, x)),
      ci = paste0("[", conf.low, ", ", conf.high, "]"),
      predictor = pred_labels
    ) |> 
    dplyr::select(predictor, estimate, ci, statistic, p)
}


# collate regtest info

tidy_regtest <- function(regtest, pb = F){
  tibble(
    estimate = ifelse(pb, regtest$b, regtest$est),
    conf.low = regtest$ci.lb,
    conf.high = regtest$ci.ub,
    p.value = regtest$pval,
    statistic = regtest$zval
  )
}

collate_regtests <- function(model, digits = 3, p_digits = 3){
  fmt <- paste0("%.", digits, "f")
  
  model |>
    tidy_regtest() |> 
    mutate(
      across(where(is.numeric), \(x) sprintf(fmt = fmt, x)),
      ci = paste0("[", conf.low, ", ", conf.high, "]")
    ) |> 
      select(estimate, ci, statistic, p.value)
}


#' Compute the full sample variance/SD from the means, SDs and sample sizes of two subgroups
#' 
#' `pooled_var()` is a helper function that computes a full sample variance estimate based on means and variances from two subgroups. This is done using the following equation
#' \deqn{s^2_\text{combined} = \frac{1}{n_1 + n_2 -1}\left((n_1-1)s_1^2 + (n_2-1)s_2^2 + \frac{n_1n_2}{n_1+n_2}(\overline{x}_1-\overline{x}_2)^2\right)}
#' This equation is essentially the same as *Result 1* in O'Neill (2014).
#' 
#' 
#' @param nx numeric value of the sample size for group 1
#' @param ny numeric value of the sample size for group 2
#' @param sdx numeric value of the standard deviation for group 1
#' @param sdy numeric value of the standard deviation for group 2
#' @param ux numeric value of the mean for group 1
#' @param uy numeric value of the mean for group 2
#' @param sd If `FALSE` returns the variance, if `TRUE` returns the standard deviation
#' 
#' @return The combined standard deviation or variance of the two groups.
#' @references 
#' O’Neill, B. (2014). Some useful moment results in sampling problems. *The American Statistician*, 68(4), 282–296. [https://doi.org/10.1080/00031305.2014.966589](https://doi.org/10.1080/00031305.2014.966589)
#' 

# function to compute full sample SD using O'Neill (2014) result 1 (variance decomposition)

pooled_var <- function(nx, ny, sdx, sdy, ux, uy, sd = F){
  n <- nx + ny
  uxy <- (nx*ux + ny*uy)/n
  pooled_var <- ((nx-1)*sdx^2 + (ny-1)*sdy^2 + (nx*ny*(ux-uy)^2)/n)/(n-1)
  if(sd){
    sqrt(pooled_var)
  } else {
    pooled_var
  }
}

#' Apply correction to Cohen's *d* to make it the unbiased Hedges \eqn{g^\ast}
#' 
#' Although when applying `escalc` from `metafor` *d* is converted to \eqn{g^\ast}, when the raw information is not available to compute *d*
#' it may be necessary to find it indirectly and convert to Hedges statistic manually. `d_to_g()` does this adjustment. 
#' 
#' @param n the sample size on which *d* is based
#' @param d the value of *d*
#' 
#' @return Hedges \eqn{g^\ast}
#' @references 
#' Hedges, L. (1981). Distribution Theory for Glass’s Estimator of Effect Size and Related Estimators. *Journal of Educational Statistics*, *6*, 107–28.
#' 


d_to_g <- function(n, d){
  (1-3/(4*(n)-9))*d
}

#' Estimate Cohen's *d* based on a correlation coefficient
#' 
#' When the raw information is not available to compute *d*, #' it may be necessary to find it indirectly. 
#' These conversions are approximations only. `d_from_r()` estimates Cohen's *d* from a correlation coefficient.
#' When `biserial = TRUE` the standard formula for converting a point biserial correlation is used, when `biserial = FALSE` the function
#' treats the measures as continuous. In this case two options need to be set based on uses Mathur & VanderWeele (2019). First, `delta`
#' is the unit change required and should be set to the mean difference between groups on the key continuous measure, *X*). second, `sx` must be set as the 
#' standard deviation of the key continuous measure, *X*). *X* is the continuous variable version of the variable on which *d* is computed for groups. 
#' For example, if *d* compares anxious vs controls, *X* is a continuous measure of anxiety.
#' 
#' @param delta the unit change required. This should be a numeric value set as the mean difference between groups on the key continuous measure, *X*.
#' @param sx   numeric value of the standard deviation of the key continuous measure, *X*
#' @param r   The correlation coeffcient to be converted
#' @param biserial When `TRUE` the standard formula for converting a point biserial correlation is used, when `FALSE` the function
#' treats the measures as continuous. In this case two options need to be set based on uses Mathur & VanderWeele (2019)
#' 
#' @return Cohen's *d*
#' @references 
#' Mathur, M. B., & VanderWeele, T. J. (2020). A simple, interpretable conversion from Pearson's correlation to Cohen's *d* for continuous exposures. *Epidemiology*, 31(2), e16–e17. [https://doi.org/10.1097/EDE.0000000000001111](https://doi.org/10.1097/EDE.0000000000001111)
#'


d_from_r <- function(delta, sx, r, biserial = T){
  if(biserial){
    2*r/(sqrt(1-r^2))
  } else {
    delta*r/(sx*sqrt(1-r^2))
  }
}

#' Estimate the sampling variance of Cohen's *d* based on a correlation coefficient
#' 
#' When the raw information is not available to compute *d*, #' it may be necessary to find it indirectly. 
#' These conversions are approximations only. [`d_from_r()`] estimates Cohen's *d* from a correlation coefficient and `var_d_from_r()` estimates the
#' associated sampling variance.
#' 
#' When `biserial = TRUE` the standard formula for converting a point biserial correlation is used, when `biserial = FALSE` the function
#' treats the measures as continuous. In this case two options need to be set based on uses Mathur & VanderWeele (2019). First, `delta`
#' is the unit change required and should be set to the mean difference between groups on the key continuous measure, *X*). second, `sx` must be set as the 
#' standard deviation of the key continuous measure, *X*). *X* is the continuous variable version of the variable on which *d* is computed for groups. 
#' For example, if *d* compares anxious vs controls, *X* is a continuous measure of anxiety.
#' 
#' @param n the total sample size on which *d* was based
#' @param r   The correlation coefficient on which *d* was based
#' @param d   The estimate of *d*
#' @param biserial When `TRUE` the standard formula for converting a point biserial correlation is used, when `FALSE` the function
#' treats the measures as continuous. In this case two options need to be set based on uses Mathur & VanderWeele (2019)
#' @param var When `TRUE` the sampling variance is returned, when `FALSE` the stancdard error is returned.
#' 
#' @return Sampling variance or standard error of Cohen's *d*
#' @references 
#' Mathur, M. B., & VanderWeele, T. J. (2020). A simple, interpretable conversion from Pearson's correlation to Cohen's *d* for continuous exposures. *Epidemiology*, 31(2), e16–e17. [https://doi.org/10.1097/EDE.0000000000001111](https://doi.org/10.1097/EDE.0000000000001111)
#'

var_d_from_r <- function(n, r, d, biserial = T, var = T){
  if(biserial){
    se <- 2/sqrt((n-1)*(1-r^2))
  } else {
    se <- abs(d)*sqrt(1/(r^2*(n-3)) + 1/(2*(n-1)))
  }
  
  if(var){
    se^2
  } else {
    se
  }
}


get_bubble_data <- function(model){
  data <- model |> metafor::regplot()
  predict <- predict(model)
  tibble::tibble(
    xi = data$xi,
    yi = data$yi,
    psize = data$psize,
    pred = predict$pred,
    pil = predict$pi.lb,
    piu = predict$pi.ub,
    cil = predict$ci.lb,
    ciu = predict$ci.ub
  )
}

#' Bubble plots from a `regtest()` object
#' 
#' `plot_bubble()` is a helper function to create a bubble plot (using `ggplot2`) based on a `regtest()` object.
#' 
#' @param model A meta-analyses object (`rma` or `rma.mv`) created using the `metafor` package
#' @param xlim  Limits for the *x*-axis
#' @param ylim  Limits for the *y*-axis
#' @param xbreaks  Breaks for the *x*-axis
#' @param ybreaks  Breaks for the *y*-axis
#' @param xlab  Text label for the *x*-axis
#' @param ylab   Text label for the *y*-axis
#' @param title  Text title for the plot
#' @param pi If `TRUE` uses a 95% prediction interval, if `FALSE` uses a 95% confidence interval
#' 
#' @return A tibble containing (for each predictor in the model) the name of the predictor (obtained from `pred_labels`), the effect size estimate, 95% confidence interval, test statistic and *p*-value.
#'


plot_bubble <- function(model, xlim = NULL, ylim = NULL, xbreaks = NULL, ybreaks = NULL, xlab = "Predictor", ylab = "Effect size", title = NULL, pi = TRUE){
  tibble <- get_bubble_data(model)
  gg <- ggplot2::ggplot(tibble, aes(x = xi, y = yi)) +
    geom_point(aes(size = psize), fill = "grey80", alpha = 0.2) +
    geom_hline(yintercept = 0, colour = "grey40") +
    geom_line(aes(x = xi, y = pred), colour = "grey10") +
    labs(x = xlab, y = ylab, title = title) +
    theme_minimal() +
    theme(legend.position = "none")
  
  if(pi){
    gg <- gg +
      geom_line(aes(x = xi, y = pil), colour = "grey10", linetype = "longdash", linewidth = 0.25) +
      geom_line(aes(x = xi, y = piu), colour = "grey10", linetype = "longdash", linewidth = 0.25)
  } else {
    gg <- gg +
      geom_line(aes(x = xi, y = cil), colour = "grey10", linetype = "longdash", linewidth = 0.25) +
      geom_line(aes(x = xi, y = ciu), colour = "grey10", linetype = "longdash", linewidth = 0.25)
  }
  
  if(!is.null(xlim)){
    if(!is.null(ylim)){
      gg <- gg + coord_cartesian(xlim = xlim, ylim = ylim)
    } else {
      gg <- gg + coord_cartesian(xlim = xlim)
    }
  } else {
    if(!is.null(ylim)){
      gg <- gg + coord_cartesian(ylim = ylim)
    }
  }
  
  try(
    if(!is.null(xbreaks)){
      gg <- gg + scale_x_continuous(breaks = seq(xlim[1], xlim[2], xbreaks))
    },
    stop("You must set xlim to use xbreaks")
  )
  
  try(
    if(!is.null(ybreaks)){
      gg <- gg + scale_y_continuous(breaks = seq(ylim[1], ylim[2], ybreaks))
    },
    stop("You must set ylim to use ybreaks")
  )
  
  gg
}

## helper functions for publication bias
# Collates models of moderate nad severe publication bias across levels of a factor and tabulates them

collate_pub_bias <- function(model, pb_mod, pb_sev, digits = 3, p_digits = 3){
  fmt <- paste0("%.", digits, "f")
  
  pb_mod <- tidy(pb_mod, conf.int = T) |> 
    select(estimate, conf.low, conf.high, p.value)
  
  pb_sev <- tidy(pb_sev, conf.int = T) |> 
    select(estimate, conf.low, conf.high, p.value)
  
  out <- pb_mod |> 
    bind_rows(pb_sev) |> 
    mutate(theta = model$b[1],
           pb = c("Moderate", "Severe")) |> 
    mutate(
      p.value = metafor::fmtp(p.value, p_digits),
      across(where(is.numeric), \(x) sprintf(fmt = fmt, x)),
      ci = paste0("[", conf.low, ", ", conf.high, "]")
    ) |> 
    select(pb, theta, estimate, ci, p.value)
  
  out
}



#' Publication bias weight models for sub-categories
#' 
#' `get_pbm()` is a helper function to fit and collate publication bias models across categories of a predictor variable. It is assumed that you will supply 
#' two vectors of values one representing moderate publication bias and the other representing severe.
#' For more information see [`metafor::selmodel()`]
#' 
#' @param tibble A tibble containing the raw data
#' @param predictor  The name of the categorical predictor for which you want individual weight models
#' @param digits  number of decimal places to print in the output
#' @param p_digits  number of decimal places for *p*-values
#' @param a numeric vector of one or more values that will be used to set the `steps` argument in `metafor::selmodel()`
#' @param moderate A vector of weights for moderate publication bias (eg. `c(1, .7, .6, .5)`) that will be used to set the `delta` argument in `metafor::selmodel()`
#' @param severe A vector of weights for severe publication bias (eg. `c(1, .5, .4, .1)`) that will be used to set the `delta` argument in `metafor::selmodel()`
#' 
#' @return A tibble containing (for each category of the predictor) a row reporting the moderate publication bias results and another reporting the severe. The table contains the 
#' original effect size estimate (\eqn{\hat{\theta}}) the estimate adjusted for bias ((\eqn{\hat{\theta}_{\text{Adj}}})), the 95% confidence interval for the adjusted value and its associated *p*-value.
#' This tibble can be printed using `kable()` in quarto.
#'



get_pbm <- function(tibble, predictor, digits = 2, p_digits = 3, a, moderate, severe){
  mas <- tibble  |>
    dplyr::arrange({{predictor}}) |> 
    dplyr::group_by({{predictor}})  |>  
    tidyr::nest()  |> 
    dplyr::mutate(
      model = purrr::map(.x = data, 
                         .f = \(es_tib) rma(yi = g, vi = v_g, data = es_tib)),
      pb_mod = purrr::map(.x = model, 
                          .f = \(m) selmodel(m, type = "stepfun", steps = a, delta = moderate)),
      pb_sev = purrr::map(.x = model, 
                          .f = \(m) selmodel(m, type = "stepfun", steps = a, delta = severe)),
      coefs = purrr::pmap(.l = list(model, pb_mod, pb_sev),
                          .f = \(model, pb_mod, pb_sev) collate_pub_bias(model, pb_mod, pb_sev, digits = digits, p_digits = p_digits)))
  
  mas
}