<div><img src="logo_ETN.png" width="100px" align="right"></div>
<div><img src="logo_share.jpg" width="100px" align="right"></div>

<div><img src="logo.png" width="200px" align="left"></div>

---
# Visualising Generalised Additive Models (GAMs). An introduction to the package gratia.
### _Rubén Rabaneda Bueno_

---


## Table of Contents

**[1. Introduction to gratia](#heading--1)**

**[2. Preparing the data](#heading--2)**

  * [2.1. Load the dataset](#heading--2-1)
  * [2.2. Fitted GAM models](#heading--2-2)
  * [2.3. Export the GAM model output to an HTML table](#heading--2-3)

**[3. Visualization of GAMs using the gratia package](#heading--3)**

  * [3.1. Inspection and visualization of model diagnostics](#heading--3-1)

      * [3.1.1. Histogram of the model residuals](#heading--3-1-1)
      * [3.1.2. Residuals versus linear predictor values](#heading--3-1-2)
      * [3.1.3. Quantile-quantile of model residuals](#heading--3-1-3)
      * [3.1.4. Fitted against observed response values](#heading--3-1-4)
      * [3.1.5. Complete model diagnostic](#heading--3-1-5)
      * [3.1.6. Rootogram of model diagnostic](#heading--3-1-6)

  * [3.2. Visualization of the main smooth effects of a GAM model](#heading--3-2)

      * [3.2.1. Estimated smooth effects](#heading--3-2-1)
      * [3.2.2. Estimated smooths effects of covariates](#heading--3-2-2)
      * [3.2.3. Random-effects smooths](#heading--3-2-3)
      * [3.2.4. Dimensional splines and tensor product smooths](#heading--3-2-4)
      * [3.2.5. Posterior smooths](#heading--3-2-5)
      * [3.2.6. Differences in smooth effects between two models](#heading--3-2-6)
      * [3.2.7. Difference smooth of factor-by smooth interactions](#heading--3-2-7)
      * [3.2.8. Penalty matrices of smooths](#heading--3-2-8)
      * [3.2.9. Derivatives of smooths](#heading--3-2-9)

**[4. Other packages for GAM visualization](#heading--4)**

**[5. Exercises](#heading--5)**

----


## <a name="heading--1"></a>1. Introduction to gratia

___gratia (Graceful 'ggplot'-Based Graphics and Other Functions for GAMs Fitted Using 'mgcv')___ is an R library developed by Gavin Simpson. As the name in parentheses suggests, it is the most _ggplot_-friendly package available for visualizing GAM models.

At the bottom of this page you will find alternative examples of graphical representations of these models, as well as some exercises that may help you become familiar with converting and formatting graphics for publication.

The best introduction to using _gratia_ to visualise the effects of GAM models is, as with all applications of R, hands-on. Let us first load the data and start with some libraries that can be useful for extracting information from GAM models in the form of tables.

## <a name="heading--2"></a>2. Preparing the data

### <a name="heading--2-1"></a>2.1. Load the dataset

__Note__ _that you can click on the icon to the left of the library load command to go directly to the library´s CRAN page._

[:books:](https://cran.r-project.org/web/packages/tidyverse/index.html)`library(tidyverse)`     # _load the readr() function_

[:books:](https://cran.r-project.org/web/packages/data.table/index.html)`library(data.table)`

[:books:](https://cran.r-project.org/web/packages/zoo/index.html)`library(zoo)`                 # _to avoid datetime format errors_

``` r
data_range <- data.table(read_csv("./data_range.csv"))
names(data_range)
```
``` r
 [1] "fishID"     "species"    "date"       "time"       "week_num"   "month_num"  "season"     "season_num" "longit_r"   "bl"
[11] "temp"
```

``` r
data_range <- as.data.frame(data_range)
data_range$time<-as.integer(data_range$time)                    # to fit the model, convert to integer
#data_range$date <- as.Date(data_range$date)
data_range$fishID<-as.factor(data_range$fishID)
data_range$species<-as.factor(data_range$species)
data_range$season<-as.factor(data_range$season)
```

``` r
head(data_range, 6)
```
| fishID 	| species 	| date 	| time 	| week_num 	| month_num 	| season 	| season_num 	| longit_r 	| bl 	| temp 	|
|---	|---	|---	|---	|---	|---	|---	|---	|---	|---	|---	|
| T449310_1 	| pikeperch 	| 27/04/2017 	| 17283 	| 5 	| 4 	| spring_I 	| 1 	| 0 	| 405 	| 8.990546 	|
| T449317_1 	| pikeperch 	| 27/04/2017 	| 17283 	| 5 	| 4 	| spring_I 	| 1 	| 536.8212 	| 415 	| 8.990546 	|
| T449202_1 	| pikeperch 	| 27/04/2017 	| 17283 	| 5 	| 4 	| spring_I 	| 1 	| 682.8973 	| 430 	| 8.990546 	|
| T449207_1 	| pike 	| 27/04/2017 	| 17283 	| 5 	| 4 	| spring_I 	| 1 	| 266.7728 	| 500 	| 8.990546 	|
| T449287_1 	| pike 	| 27/04/2017 	| 17283 	| 5 	| 4 	| spring_I 	| 1 	| 223.4977 	| 595 	| 8.990546 	|
| T449313_1 	| pikeperch 	| 27/04/2017 	| 17283 	| 5 	| 4 	| spring_I 	| 1 	| 0 	| 495 	| 8.990546 	|

### <a name="heading--2-2"></a>2.2. Fitted GAM models

[:books:](https://cran.r-project.org/web/packages/mgcv/index.html)`library(mgcv)`

The models used here are GAM with mixed effects and an AR autocorrelation structure (2). All but one of the models were fitted using the _bam()_ function and include species as a parametric term and as a smooth factor variable.

``` r
mod <- bam(longit_r ~ species +
                      s(time, by = species, bs = "cc") +
                      s(fishID, bs = 're'),        # same as s(species, fishID, bs = 're')
                      data = data_range,
                      family = 'gaussian', correlation = corARMA(form = ~ 1|time, p = 2))
```

Other models used for specific purposes are shown after clicking on the black arrow:

<details><summary>Model <strong>"mod2"</strong> (t.p.r.s. smooths)</summary><blockquote>

``` r
mod2 <- bam(longit_r ~ species +
                       s(time, by = species, bs = "ts") +
                       s(fishID, bs = 're'),        # same as s(species, fishID, bs = 're')
                       data = data_range,
                       family = 'gaussian', correlation = corARMA(form = ~ 1|time, p = 2))
```

  </blockquote></details>
</blockquote></details>

<details><summary>Model <strong>"mod_fs"</strong> (factor smooth interaction)</summary><blockquote>

``` r
mod_fs <- bam(longit_r ~ species +
                         s(temp) + s(bl)+
                         s(temp, fishID, bs = "fs", k=10, m = 1) +
                         s(bl, fishID, bs = "fs", k=10, m = 1) +
                         s(time, by = species, bs = "cc") +
                         s(time, fishID, bs = "fs", k=10, m=1),
                         data = data_range,
                         family = 'gaussian', correlation = corARMA(form = ~ 1|time, p = 2))
```

  </blockquote></details>
</blockquote></details>

<details><summary>Model <strong>"mod_tps"</strong> (interacting t.p.r.s.)</summary><blockquote>

``` r
mod_tps <- bam(longit_r ~ species +
                          s(temp, week_num, by = species) +
                          s(fishID, bs = 're'),
                          data = data_range,
                          family = 'gaussian', correlation = corARMA(form = ~ 1|time, p = 2))
```

<details><summary>Model <strong>"mod_te"</strong> (tensor product smooths)</summary><blockquote>

``` r
mod_te <- bam(longit_r ~  species +
                          s(time, by = species, bs="cc") +
                          te(bl, temp, by = species),
                          data = data_range,
                          family = 'gaussian', correlation = corARMA(form = ~ 1|time, p = 2))
```

  </blockquote></details>
</blockquote></details>

<details><summary>Model <strong>"mod_pike"</strong> (gamm/pike)</summary><blockquote>

``` r
mod_pike <- gamm(longit_r ~ s(time, bs = "cc", k = 12),
                            random = list(fishID = ~1),
                            family = 'gaussian',
                            data = data_range %>%
                            filter(species == "pike") %>%
                            mutate(fishID = as_factor(fishID)),
                            correlation = corARMA(form = ~ 1|time, p = 2))
```

  </blockquote></details>
</blockquote></details>


### <a name="heading--2-3"></a>2.3. Export the GAM model output to an HTML table

- Using the _tab_model()_ function from the _sjPlot_ package.

[:books:](https://cran.r-project.org/web/packages/sjPlot/index.html)`library(sjPlot)`

``` r
sjPlot::tab_model(mod,
                      title='GAM model 1 (mod)',
#                     pred.labels = c('Intercept',
#                                     'species[Pikeperch]',
#                                     'species[Wels]',
#                                     'smooth(Time[Pike])',
#                                     'smooth(Time[Pikeperch])',
#                                     'smooth(Time[Wels])',
#                                     'fishID'),
                     string.pred = "Parameter",
                     string.p = "P-Value",       # Use p.style = "stars" for significance levels as stars instead
                     show.ci=FALSE
)
```

- Using the _gamtabs()_ function from the _itsadug_ package.

[:books:](https://cran.r-project.org/web/packages/itsadug/index.html)`library(itsadug)`

``` r
itsadug::gamtabs(mod, caption="", comment=FALSE, type='html')            # use $gam for gamm4 objet
```

<table border=1>
<caption align="bottom">  </caption>
  <tr> <td> A. parametric coefficients </td> <td align="right"> Estimate </td> <td align="right"> Std. Error </td> <td align="right"> t-value </td> <td align="right"> p-value </td> </tr>
  <tr> <td> (Intercept) </td> <td align="right"> 764.2741 </td> <td align="right"> 144.6854 </td> <td align="right"> 5.2823 </td> <td align="right"> &lt; 0.0001 </td> </tr>
  <tr> <td> speciespikeperch </td> <td align="right"> 328.0665 </td> <td align="right"> 218.5919 </td> <td align="right"> 1.5008 </td> <td align="right"> 0.1335 </td> </tr>
  <tr> <td> specieswels </td> <td align="right"> 432.0166 </td> <td align="right"> 185.3349 </td> <td align="right"> 2.3310 </td> <td align="right"> 0.0198 </td> </tr>
   <tr> <td> B. smooth terms </td> <td align="right"> edf </td> <td align="right"> Ref.df </td> <td align="right"> F-value </td> <td align="right"> p-value </td> </tr>
  <tr> <td> s(time):speciespike </td> <td align="right"> 1.0376 </td> <td align="right"> 8.0000 </td> <td align="right"> 0.2274 </td> <td align="right"> 0.1736 </td> </tr>
  <tr> <td> s(time):speciespikeperch </td> <td align="right"> 5.1923 </td> <td align="right"> 8.0000 </td> <td align="right"> 5.4342 </td> <td align="right"> &lt; 0.0001 </td> </tr>
  <tr> <td> s(time):specieswels </td> <td align="right"> 5.9984 </td> <td align="right"> 8.0000 </td> <td align="right"> 16.7412 </td> <td align="right"> &lt; 0.0001 </td> </tr>
  <tr> <td> s(fishID) </td> <td align="right"> 25.6476 </td> <td align="right"> 27.0000 </td> <td align="right"> 18.7155 </td> <td align="right"> &lt; 0.0001 </td> </tr>
   <a name=tab.gam></a>
</table>

__Tip__: For fitting with gamm use the gam component of the fitted model

## <a name="heading--3"></a>3. Visualization of GAMs with the _gratia_ package

[:books:](https://cran.r-project.org/web/packages/gratia/index.html)`library(gratia)`

In the next sections we will use several functions from this package, but also some additional functions that are deprecated (including `derivSimulCI()`, `plot.derivSimulCI()`, `signifD()`) and that we download from [G. Simpson´s GitHub repository](https://gist.github.com/gavinsimpson/ca18c9c789ef5237dbc6).

Download additional functions not currently included in the updated version of _gratia_.

[:books:](https://cran.r-project.org/web/packages/curl/index.html)`library(curl)`

``` r
tmpf <- tempfile()
curl_download("https://gist.githubusercontent.com/gavinsimpson/ca18c9c789ef5237dbc6/raw/295fc5cf7366c831ab166efaee42093a80622fa8/derivSimulCI.R", tmpf)
source(tmpf)
```
<details><summary><strong>Show functions code</strong></summary><blockquote>

``` r
`derivSimulCI` <- function(mod, n = 200, eps = 1e-7, newdata, term,
                           samples = 10000) {
    stopifnot(require("MASS"))
    if(inherits(mod, "gamm"))
        mod <- mod$gam
    m.terms <- attr(terms(mod), "term.labels")
    if(missing(newdata)) {
        newD <- sapply(model.frame(mod)[, m.terms, drop = FALSE],
                       function(x) seq(min(x), max(x) - (2*eps), length = n))
        names(newD) <- m.terms
    } else {
        newD <- newdata
    }
    newDF <- data.frame(newD) ## needs to be a data frame for predict
    X0 <- predict(mod, newDF, type = "lpmatrix")
    newDF <- newDF + eps
    X1 <- predict(mod, newDF, type = "lpmatrix")
    Xp <- (X1 - X0) / eps
    Xp.r <- NROW(Xp)
    Xp.c <- NCOL(Xp)
    ## dims of bs
    bs.dims <- sapply(mod$smooth, "[[", "bs.dim") - 1
    ## number of smooth terms
    t.labs <- attr(mod$terms, "term.labels")
    ## match the term with the the terms in the model
    if(!missing(term)) {
        want <- grep(term, t.labs)
        if(!identical(length(want), length(term)))
            stop("One or more 'term's not found in model!")
        t.labs <- t.labs[want]
    }
    nt <- length(t.labs)
    ## list to hold the derivatives
    lD <- vector(mode = "list", length = nt)
    names(lD) <- t.labs
    ## sample draws from the posterior distribution of model coefficients
    Rbeta <- t(mvrnorm(n = samples, coef(mod), vcov(mod)))
    ## loop over the terms
    for(i in seq_len(nt)) {
        want <- grep(t.labs[i], colnames(X1))
        lD[[i]] <- list(deriv = Xp[, want] %*% coef(mod)[want],
                        simulations = Xp[, want] %*% Rbeta[want, ])
    }
    class(lD) <- "derivSimulCI"
    lD$gamModel <- mod
    lD$eps <- eps
    lD$eval <- newD - eps
    lD ##return
}

plot.derivSimulCI <- function(x, alpha = 0.05, polygon = TRUE,
                              sizer = FALSE, term,
                              eval = 0, lwd = 3,
                              col = "lightgrey", border = col,
                              ylab, xlab, main, ...) {
    l <- length(x) - 3
    ## get terms and check specified (if any) are in model
    term.labs <- names(x[seq_len(l)])
    if(missing(term)) {
        term <- term.labs
    } else {
        term <- term.labs[match(term, term.labs)]
    }
    if(any(miss <- is.na(term)))
        stop(paste("'term'", term[miss], "not a valid model term."))
    if(all(miss))
        stop("All terms in 'term' not found in model.")
    l <- sum(!miss)
    nplt <- n2mfrow(l)
    if(missing(ylab))
        ylab <- expression(italic(hat(f)*"'"*(x)))
    if(missing(xlab)) {
        xlab <- attr(terms(x$gamModel), "term.labels")
        names(xlab) <- xlab
    }
    if (missing(main)) {
        main <- term
        names(main) <- term
    }
    ## compute confidence interval
    ciFUN <- function(x, alpha) {
        ahalf <- alpha / 2
        apply(x$simulations, 1, quantile, probs = c(ahalf, 1 - ahalf))
    }
    CI <- lapply(x[seq_len(l)], ciFUN, alpha = alpha)
    ## plots
    layout(matrix(seq_len(l), nrow = nplt[1], ncol = nplt[2]))
    on.exit(layout(1))
    for(i in term) {
        lwr <- CI[[i]][1,]
        upr <- CI[[i]][2,]
        ylim <- range(upr, lwr)
        plot(x$eval[,i], x[[i]]$deriv, type = "n",
             ylim = ylim, ylab = ylab, xlab = xlab[i], main = main[i], ...)
        if(isTRUE(polygon)) {
            polygon(c(x$eval[,i], rev(x$eval[,i])),
                    c(upr, rev(lwr)), col = col, border = border)
        } else {
            lines(x$eval[,i], upr, lty = "dashed")
            lines(x$eval[,i], lwr, lty = "dashed")
        }
        abline(h = 0, ...)
        if(isTRUE(sizer)) {
            lines(x$eval[,i], x[[i]]$deriv, lwd = 1)
            S <- signifD(x[[i]]$deriv, x[[i]]$deriv, upr, lwr,
                         eval = eval)
            lines(x$eval[,i], S$incr, lwd = lwd, col = "blue")
            lines(x$eval[,i], S$decr, lwd = lwd, col = "red")
        } else {
            lines(x$eval[,i], x[[i]]$deriv, lwd = 2)
        }
    }
    invisible(x)
}

signifD <- function(x, d, upper, lower, eval = 0) {
    miss <- upper > eval & lower < eval
    incr <- decr <- x
    want <- d > eval
    incr[!want | miss] <- NA
    want <- d < eval
    decr[!want | miss] <- NA
    list(incr = incr, decr = decr)
}
```

``` r
tmpf <- tempfile()
curl_download("https://gist.githubusercontent.com/gavinsimpson/d23ae67e653d5bfff652/raw/25fd719c3ab699e48927e286934045622d33b3bf/simulate.gamm.R", tmpf)
source(tmpf)
```

  </blockquote></details>
</blockquote></details>


### <a name="heading--3-1">3.1. Inspection and visualization of model diagnostics

### <a name="heading--3-1-1">3.1.1. Histogram of the model residuals

Plot the histogram of the model residuals with the function  `gratia::residuals_hist_plot()`

``` r
plt_res_hist <- gratia::residuals_hist_plot(mod, type = "pearson",                               # other types of residuals include "deviance" and "response"
                                            title = "Pearson residuals of the fitted GAM",
                                            subtitle = "mod 1",
                                            caption = NULL)
```

One of the advantages of the package is that its implemented functions are based on the layer system of the _ggplot2_ package and are therefore easily customizable. For example, use a custom theme for the size and colour of the fonts of axes and titles.

[:books:](https://cran.r-project.org/web/packages/ggplot2/index.html)`library(ggplot2)`

Define new gpplot themes and colour palette.

``` r
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")    # define a palette color

th1 <- theme_bw() +
        theme(
              plot.title = element_text(color="red", size=16, face="bold"),
              plot.subtitle = element_text(color="black", size=12, face="italic"),
              axis.title.x = element_text(size = 14),
              axis.text.x = element_text(size = 14),
              axis.title.y = element_text(size = 14),
              axis.text.y = element_text(size = 14),
              axis.text = element_text(size = 14))

th2 <- theme_bw() &
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              plot.title = element_text(color="black", size=16, face="bold"),
              plot.subtitle = element_text(color="black", size=16, face="italic"),
              axis.title.x = element_text(size = 16),
              axis.text.x = element_text(size = 16),
              axis.title.y = element_text(size = 16),
              axis.text.y = element_text(size = 16),
              axis.text = element_text(size = 16))
```

Draw the histogram of the residuals.
``` r
plt_res_hist + th1
```

<p align="center"><img src="/Plots/Fig.res_hist.png?raw=true"></p>

### <a name="heading--3-1-2">3.1.2. Residuals versus linear predictor values

Plot residuals versus linear predictor values with the function `gratia::residuals_linpred_plot()`

``` r
plt_res_pred <- gratia::residuals_linpred_plot(mod, type = "deviance",
                                                    point_col = cbPalette[1],
                                                    point_alpha = 0.5,
                                                    line_col = cbPalette[4])

plt_res_pred + th1
```

<p align="center"><img src="/Plots/Fig.res_pred.png?raw=true"></p>

### <a name="heading--3-1-3">3.1.3. Quantile-quantile of model residuals

To draw Q-Q plots of model residuals use the function `qq_plot()`

The _method_ parameter allows us to specify how the theoretical quantiles are generated. The default method "uniform" uses direct randomization, while `simulate`performs simulations using the model and `uniform` represents the normality assumption.

Alternatively, the `gratia::worm_plot()` function can be used to check the model residuals, specifying the same method parameters as for the `qq_plot()` function.

``` r
plt_qq <- gratia::qq_plot(mod, method = "uniform",
                               point_col = cbPalette[4],
                               point_alpha = 0.4)

plt_qq + th1
```

<p align="center"><img src="/Plots/Fig.plt_qq.png?raw=true"></p>

### <a name="heading--3-1-4">3.1.4. Fitted against observed response values

Plot fitted against observed response values with the function `observed_fitted_plot()`

``` r
gratia::observed_fitted_plot(mod, point_col = cbPalette[1],
                                  point_alpha = 0.5) +
                                  th1
```

<p align="center"><img src="/Plots/Fig.obs_vs_fit.png?raw=true"></p>

### <a name="heading--3-1-5">3.1.5. Complete model diagnostic

Instead of using each of the above functions individually, we can use the function `appraise()` to draw a comprehensive model diagnosis.

``` r
gratia::appraise(mod, type = c("pearson"),
                      method = c("uniform"),
                      point_col = cbPalette[1],
                      point_alpha = 0.5,
                      line_col = cbPalette[4]) &
                      th2
```

<p align="center"><img src="/Plots/Fig.complete_diag.png?raw=true"></p>

### <a name="heading--3-1-6">3.1.6. Rootogram of model diagnostic

Draw a rootogram model diagnostic with the functions `gratia::rootogram()`. The rootogram is used to assess the godness-of-fit of the model by comparing observed to expected frequencies of response values.

``` r
root_mod <- gratia::rootogram(mod)

gratia::draw(root_mod, type = "standing") & th1              # the other types of rootgram are "hanging" and "suspended"
```

<p align="center"><img src="/Plots/Fig.rootg.png?raw=true"></p>

### <a name="heading--3-2">3.2. Visualization of the main smooth effects of a GAM model

### <a name="heading--3-2-1">3.2.1. Estimated smooth effects

The main function for drawing GAM models is `gratia::draw()`, which is based on the function `ggplot2::ggplot()`.

The output of most gratia functions is a tibble table which, which we can manipulate either for further plot adjustments with _ggplot()_ or for data analysis.

``` r
gratia::draw(mod, select = NULL,        # plot all smooth estimates
                  residuals = FALSE,    # do not add partial residuals
                  scales = "free",      # all smooths in own y-axis scale. If not, "fixed"
#                 unconditional = TRUE, # if TRUE, corrected Bayesian covar. matrix is used
#                 n = NULL,             # num. points of covar. range on to evaluate smoothing
#                 ci_col = "black",     # set color for CIs
#                 smooth_col = "white", # set color for smooth or difference curve
#                 ci_alpha = 0.3        # adjust transparency for CIs
                  ) & th2               # add a ggplot theme with '&'
```
<p align="center"><img src="/Plots/Fig.sm_eff.png?raw=true"></p>

- The plot shows the quantiles of the smoothing function for random effects (i.e., the fish-specific intercepts fitted with `bs='re'`) and time trends of different shape and wiggliness for each species with the factor-by smooth (`by`).

- With the parameter `select` we can choose the smooth term of the model to be drawn and with `unconditional` we choose whether the uncertainty due to smoothness selection is included in the CIs.

### <a name="heading--3-2-2">3.2.2. Evaluation of the estimated smooths of covariates

The parameter _`n`_ is used to evaluate the smooths of covariates on a grid of n points uniformly distributed over their range.

It is implemented in the function `smooth_estimates()`, whose tibble output contains the values of the smoothing for the desired range of the grid n points over one or all covariates associated with the smooth.

- With the parameter `smooth` we can choose the smooth term of the model to evaluate. For example, to evaluate the smooth of time of the model __mod__ at 20 points (albeit suboptimal!):

``` r
sm_eval <- gratia::smooth_estimates(mod, smooth = "s(time)",               # choose a smooth term to evaluate. If NULL, evaluate all smooth terms.
                                         partial_match = TRUE,             # to match the smoothings
                                         n = 20)

sm_eval$time <- as.Date(sm_eval$time)                                      # convert 'time' to Date format

gratia::draw(sm_eval) & th2                                                # plot the smooths of time for each level of the 'by' variable
```
<p align="center"><img src="/Plots/Fig.sm_eval.png?raw=true"></p>

### <a name="heading--3-2-3">3.2.3. Random-effects smooths

If the model contains a random factor smooth interaction at the group level (specified as `bs='fs'`), the function 'draw()' draws the individual slopes with the same wiggle.

The __"mod_fs"__ model, unlike the __"mod"__ model, has a random smooth interaction so that instead of individual specific random effects (i.e., fish intercepts), a different slope for each fish is displayed over the range of individual covariates in addition to the smooth trends of the main effects for each species. In this case, we choose only the random smoothing factors of the covariates 'bl' and 'temp' and their smooth main effects.

``` r
gratia::draw(mod_fs, select = c(1,2,3,4),
                     guides="collect") &      # places ranef legend to the right ("keep", default, for a mid placement)
                     th2
```

<p align="center"><img src="/Plots/Fig.sm_fs.png?raw=true"></p>

### <a name="heading--3-2-4">3.2.4. Dimensional splines and tensor product smooths

If the model contains smooth surfaces due to an interaction between covariates from multidimensional thin plate splines (e.g. t.p.r.s. `s(x1, x2)`), `draw()` generates 2D views of the variables in the interaction.

``` r
gratia::draw(mod_tps, select = c(1,2,3),                                                # select only the interaction surfaces
             scales = "free",
             n_contour = 15,  n.grid=100,                                               # set the number of contour bins
#            contour_col = "white",
             continuous_fill = scale_fill_distiller(palette= "Spectral", type = "div"), # change the colour scheme (NULL is the default palette "RdBu")
             ncol = 3) &
             th2
```

<p align="center"><img src="/Plots/Fig.sm_tpr_2d.png?raw=true"></p>

If the model contains a tensor product of the smooth interaction (e.g. `te(x1, x2)`, `ti(x1, x2)`, `t2(x1, x2)`), we can also easily visualise the interacting splines with the `draw()` function. For the model __"mod_te"__, for example, we draw only the smoothed 2D surfaces with:

``` r
gratia::draw(mod_te, select = c(4,5,6),
                     scales = "free",
                     n_contour = 15,
                     ncol = 3) &
                     th2
```

<p align="center"><img src="/Plots/Fig.sm_te.png?raw=true"></p>

### <a name="heading--3-2-5">3.2.5. Posterior smooths

The function `gratia::smooth_samples()`> creates a tibble table containing draws of the posterior distributions of each smooth function in the model. The posterior smooths can then be drawn with `draw()`.

``` r
sm_post_mod <- gratia::smooth_samples(mod, n = 50,                     # computes 50 posterior draws
                                           seed = 1000)                # for reproducibility
head(sm_post_mod, 5)
```
``` r
# A tibble: 5 x 8
  smooth  term                by_variable species .x1          row  draw value
  <chr>   <chr>               <chr>       <fct>   <date>     <int> <int> <dbl>
1 s(time) s(time):speciespike species     pike    2017-04-27     1     1  36.6
2 s(time) s(time):speciespike species     pike    2017-04-27     2     1  38.3
3 s(time) s(time):speciespike species     pike    2017-04-27     3     1  40.0
4 s(time) s(time):speciespike species     pike    2017-04-28     4     1  41.6
5 s(time) s(time):speciespike species     pike    2017-04-28     5     1  43.2
```
``` r
sm_post_mod$.x1 <- as.Date(sm_post_mod$.x1)             # convert '.x1' to Date  format
```
``` r
gratia::draw(sm_post_mod, n_samples = 15,               # 'n_samples' sets the number of randomly sampled draws to diaplay
                          alpha = 0.5,
                          col = cbPalette[6],
                          nrow=1) &
                          th2
```

<p align="center"><img src="/Plots/Fig.sm_post.png?raw=true"></p>

### <a name="heading--3-2-6">3.2.6. Differences in smooth effects between two models

The function `gratia::compare_smooths()` compares smooths from different models, i.e. with different spline functions fitted to the smooth terms, whose product is a tibble table with the comparisons.

For example, let us compare the models __"mod"__, whose smooth terms were fitted with 'cc' splines, with __"mod2"__, which were fitted with 'tp' splines.

``` r
sms_comp <- gratia::compare_smooths(mod, mod2)

head(sms_comp, 5)
```
```
# A tibble: 5 x 5
  model smooth                   type          by      data
  <chr> <chr>                    <chr>         <chr>   <list>
1 mod   s(fishID)                Random effect NA      <tibble [30 × 3]>
2 mod2  s(fishID)                Random effect NA      <tibble [30 × 3]>
3 mod   s(time):speciespike      Cyclic CRS    species <tibble [100 × 4]>
4 mod2  s(time):speciespike      TPRS (shrink) species <tibble [100 × 4]>
5 mod   s(time):speciespikeperch Cyclic CRS    species <tibble [100 × 4]>
```
We can manipulate the output of difference smooth by extending the x-axis titles with a simple `for loop` to convert 'fishID' to numeric format (to avoid overlap) and convert 'time' from integer to Date format for better display.

``` r
t.dat <- list(3,4,5,6,7,8)
id.dat <- list(1,2)
for(i in t.dat){
    sms_comp$data[[i]]$time <- as.Date(sms_comp$data[[i]]$time)
  for(i in id.dat)
      {sms_comp$data[[i]]$fishID <- as.numeric(sms_comp$data[[i]]$fishID)
      }}
```
<p align="center"><img src="/Plots/Fig.sms_comp.png?raw=true"></p>

### <a name="heading--3-2-7">3.2.7. Difference smooth of factor-by smooth interactions

The function `gratia::difference_smooths()` can be used to compare the smoothing curves of the factor-by smooth variable (i.e. `by`) of a model. This function creates a tibble table with the difference smoothing between each pair of levels of the grouping variables.

If we want to compare the smoothing curves between the levels of _'species'_, we first use `gratia::difference_smooths()` to generate the difference smoothing pairs between the factor levels and then plot the differences with `draw()`>.

As with the smoothing covariates, we can evaluate the dfifference smooth for __`n`__ points of the covariates. Here, the parameter __`smooth`__ is used to select the smooth term of the model to be evaluated.

``` r
sms_diffs <- gratia::difference_smooths(mod, smooth = "s(time)", # diff. smooth for time
                                             ci_level = 0.975)   # credible interval set to 97.5%

sms_diffs$time <- as.Date(sms_diffs$time)                        # convert 'time' to Date format

head(sms_diffs, 5)
```
```
# A tibble: 5 x 9
  smooth  by      level_1 level_2    diff    se   lower upper time
  <chr>   <chr>   <chr>   <chr>     <dbl> <dbl>   <dbl> <dbl> <date>
1 s(time) species pike    pikeperch  213.  89.1  13.9    413. 2017-04-27
2 s(time) species pike    pikeperch  207.  88.4   8.33   405. 2017-04-27
3 s(time) species pike    pikeperch  196.  86.9   0.978  391. 2017-04-28
4 s(time) species pike    pikeperch  182.  85.2  -8.69   373. 2017-04-29
5 s(time) species pike    pikeperch  167.  83.8 -21.3    354. 2017-04-30
```

``` r
gratia::draw(sms_diffs, ncol = 3) & th2
```
<p align="center"><img src="/Plots/Fig.sms_diffs.png?raw=true"></p>

- We see the smooth difference between each pair of levels of the factorial smooth variable _'species'_, but `difference_smooths()`does not show regions where these differences between species are significant. A way to plot these significant periods within the smooth differences, but one that is less ggplot friendly, is to use the <span style="color:purple">`plot_diff()`</span> function from the _itsadug_ package (see alternative method below). However, We can still display these regions with <span style="color:purple">`difference_smooths()`</span> by extracting the estimated differences of the pairwise comparisons and the credible interval information from the tibble table created with the function and then calculating the periods of significant differences with the function <span style="color:purple">`signifD()`</span>.

### <a name="heading--3-2-8">3.2.8. Penalty matrices of smooths

The `gratia::penalty()` is used to extract penalty matrices of the smooths and tidy them in a tibble table, which can then be displayed as a heatmap with with `draw()`.

We can select the specific smooth for which display the penalty matrix witht the parameter `smooth` (e.g., for pike, use smooth = `'s(time):speciespike'`)

``` r
penal_mod <- gratia::penalty(mod)

gratia::draw(penal_mod, smooth = NULL) & th2
```

<p align="center"><img src="/Plots/Fig.penal_mx.png?raw=true"></p>

### <a name="heading--3-2-9">3.2.9. Derivatives of smooths

The function `derivatives()` contains a simple method for computing the first derivatives of estimated smooths from a GAM model using finite differences, where `type` specifies a `central`, `forward` or `backward` approximation.

The result is a tibble table with the estimated derivatives of all or the selected smooth and the upper and lower bounds of the specified pointwise confidence interval (CI) or simultaneous CI calculated using posterior Bayesian simulations. The derivatives can then be drawn using `draw()`.

Using the first derivative, we can evaluate the GAM spline by identifying the periods with statistically significant changes in the nonlinear trend (i.e.,the first derivatives are significantly different from zero).

If the models include factorial smoothing ('by'), you must specify the smooth on which to calculate the derivatives with the `term` parameter and set `partial_match` to `TRUE`.

Let us examine how the rate of change of the longitudinal motion behaves over time (i.e., the GAM slope) and identify the periods with significant upward or downward changes.

``` r
sms_deriv <- gratia::derivatives(mod, type = "central",     # select central finite differences
                                      term = "s(time)",     # select smooth to compute fd
                                      partial_match = TRUE) # required if GAM  incl. a `by` variable

sms_deriv$data <- as.Date(sms_deriv$data)                   # convert 'time' to Date  format

head(sms_deriv, 5)
```
```
# A tibble: 5 x 8
  smooth              var   data       derivative    se  crit lower upper
  <chr>               <chr> <date>          <dbl> <dbl> <dbl> <dbl> <dbl>
1 s(time):speciespike time  2017-04-27       1.33  2.77  1.96 -4.09  6.76
2 s(time):speciespike time  2017-04-27       1.27  2.77  1.96 -4.16  6.71
3 s(time):speciespike time  2017-04-27       1.21  2.79  1.96 -4.26  6.68
4 s(time):speciespike time  2017-04-28       1.14  2.82  1.96 -4.38  6.67
5 s(time):speciespike time  2017-04-28       1.07  2.85  1.96 -4.52  6.65
```

``` r
gratia::draw(sms_deriv, ncol = 3) & th2
```

<p align="center"><img src="/Plots/Fig.sms_deriv.png?raw=true"></p>

- Derivatives of the temporal smoothings for the three species, with the shaded area representing the 95% simultaneous CI.

Displaying periods with significant changes in fitted slope is currently not possible with <span style="color:purple">`derivatives()`</span>, but we can apply the deprecated function <span style="color:purple">`derivSimulCI()`</span> with the parameter 'sizer' set to 'TRUE' to models fitted with <span style="color:purple">`gamm()``</span> without a parametric or a factor-by-smooth variable. The first derivatives are calculated with simultaneous intervals.

For example, let us calculate the first derivatives of the pike time trend using the __mod_pike__ model.

``` r
library(curl)                      # download functions not currently incl. gratia

tmpf <- tempfile()
curl_download("https://github.com/AngeVar/GLAHD/blob/master/R/fitGAM/derivSimulCI.R", tmpf)
source(tmpf)
```
``` r
sms_deriv_pike <- derivSimulCI(mod_pike, samples = 1000,  # run simulations
                                         n = 200)
```

- The calculated derivatives are stored in a list object _derivSimulCI_, which contains the simulations used in the calculation of the simultaneous intervals.


``` r
str(sms_deriv_pike, 1)
```
```
List of 4
 $ time    :List of 2
 $ gamModel:List of 31
  ..- attr(*, "class")= chr "gam"
 $ eps     : num 1e-07
 $ eval    : num [1:200, 1] 17283 17283 17284 17284 17285 ...
  ..- attr(*, "dimnames")=List of 2
 - attr(*, "class")= chr "derivSimulCI"
 ```

With the function `plot.derivSimulCI()` we can plot and highlight the significant periods on the curve.

``` r
plot.derivSimulCI(sms_deriv_pike, sizer = TRUE)     # 'sizer' = FALSE, hides signf. periods
```
``` r
knitr::include_graphics('Fig.sms_deriv_pike.png')
```

The slope of the derivatives differs from the previous plot because the variable "species" was not included in the model as a parametric factor or "by" variable.

Periods with significant rates of change are shown in red for downward trends and in blue for upward trends (i.e., periods where the CI does not include 0).

We can adjust the previous representation in a similar way as we did for the smooth difference plot. First, we need to extract from the created list _derivSimulCI_ the first derivatives in `deriv`, the simultaneous intervals calculated using simulations in `simulations`, and the significance periods in _`incr`_ and _`decr`_.

``` r
CI <- apply(sms_deriv_pike[[1]]$simulations, 1, quantile, probs = c(0.025, 0.975))
sigD <- signifD(sms_deriv_pike[["time"]]$deriv, sms_deriv_pike[["time"]]$deriv, CI[2, ], CI[1, ],
                eval = 0)
```
``` r
sms_deriv_pike <- with(sms_deriv_pike,
                   data.frame(derivative = sms_deriv_pike[["time"]]$deriv[, 1],
                              CI_upper = (apply(sms_deriv_pike[[1]]$simulations, 1,
                                                quantile, probs = c(0.025, 0.975)))[2, ], # upper CI
                              CI_lower = (apply(sms_deriv_pike[[1]]$simulations, 1,
                                                quantile, probs = c(0.025, 0.975)))[1, ], # upper CI
                              time <- as.Date(sms_deriv_pike[[4]]),
                              Signif.incr = (signifD(sms_deriv_pike[["time"]]$deriv,
                                                     sms_deriv_pike[["time"]]$deriv,
                                                     CI[2, ], CI[1, ], eval = 0))$incr,
                              Signif.decr = (signifD(sms_deriv_pike[["time"]]$deriv,
                                                     sms_deriv_pike[["time"]]$deriv,
                                                     CI[2, ], CI[1, ], eval = 0))$decr
                                  )
                       )
colnames(sms_deriv_pike)[4] <- "time"                                                     # rename time column
```

Plot the first derivatives using the _ggplot_ package and show the significant increasing and decreasing rates of change in a different color.

``` r
ggplot(sms_deriv_pike,
       aes(x = time, y = derivative)) +
       scale_x_date(breaks = "1 month", minor_breaks = "1 day", date_labels = "%b %y") +
       geom_ribbon(aes(ymax = CI_upper, ymin = CI_lower), alpha = 0.3, fill = cbPalette[4]) +
       geom_line(col = cbPalette[4],size = 1) +
       geom_line(aes(y = Signif.incr), size = 2,col = "white") +
       geom_line(aes(y = Signif.decr), size = 2,col = "black") +
       geom_hline(data = sms_deriv_pike, aes(yintercept = 0),color="black", linetype="dashed") +
       ggtitle("Rate of change of pike longitudinal range") +
       ylab(expression(italic(hat(f) * "'") * ("estimated trend"))) +
       xlab("Month") +
       th2
```

<p align="center"><img src="/Plots/Fig.sms_deriv_pike_ggpl.png?raw=true"></p>

## <a name="heading--4">4. Other packages for GAM visualization

In addition to the _gratia_ package, there are other R packages that allow you to display various types of visualisations of GAM models, some of which are more and some of which are less ggplot-friendly.

## <a name="heading--4">4. Other packages for GAM visualization: mgcv

[:books:](https://cran.r-project.org/web/packages/mgcv/index.html)`library(mgcv)`

The function `plot.gam()` is based on R graphs, unlike `draw()`. But as with the latter, we can select the smooth term of the model with the _`select`_ parameter.

``` r
par(mfrow = c(2,2))

mgcv::plot.gam(mod, Select = NULL,
                    shade = TRUE,
                    scale = 0,
                    seWithMean = TRUE)
```
<details><summary>Show plot</summary><blockquote>

<p align="center"><img src="/Plots/Fig.sm_eff_2.png?raw=true"></p>

   </blockquote></details>
</blockquote></details>

## <a name="heading--4">4. Other packages for GAM visualization: mgcViz

[:books:](https://cran.r-project.org/web/packages/mgcViz/index.html)`library(mgcViz)`

The function `plot.gamViz()` is equivalent to `plot.gam()`, but based on the _ggplot2_ layer system and is therefore easily customizable.

To start plotting we first need to create a _gamViz_ object with the `getViz()` and then plot it with `plot.gamViz()`.

``` r
viz_mod <- mgcViz::getViz(mod)
plt_mod <- mgcViz::plot.gamViz(viz_mod, select = NULL) +              # select model terms
                                        geom_hline(yintercept = 0)    # draw horizontal line
print(plt_mod, pages = 1)
```
<details><summary>Show plot</summary><blockquote>

<p align="center"><img src="/Plots/Fig.sm_eff_3.png?raw=true"></p>

   </blockquote></details>
</blockquote></details>

_ggplot_ customization themes are usable with _gamViz_ objects.

``` r
plt_mod <- mgcViz::plot.gamViz(viz_mod, allTerms = TRUE) +        # allTerms=T displays all terms
                               l_points(col = 1, size = 0.5) +    # for scaled trend lines
                               l_fitPoints(col = 2, size=1) +     # view mean fit points
                               l_fitLine(linetype = 3, col = 4) + # line of the fitted curve
                               l_fitContour() +
                               l_ciBar() +                        # group means error bars
#                              l_ciPoly() +                       # fills shaded CIs
#                              l_dens(type = "cond")              # conditional density of p. res.
#                              l_ciLine(colour = 2) +             # draw line on CIs bands
                               labs(title = NULL) +               # labeles
                               th2
print(plt_mod, pages = 1)
```

<details><summary>Show plot</summary><blockquote>

<p align="center"><img src="/Plots/Fig.sm_eff_4.png?raw=true"></p>

   </blockquote></details>
</blockquote></details>


`plot.gamViz()()` can also be used to select the smooth term of the grouping variable with `mgcViz::sm()` and then you can combine plots with `mgcViz::gridPrint()`.

To combine the smooth trend plots of wels from the __"mod"__ and __"mod2"__ models, first select the corresponding smooth terms.

``` r
viz_mod2 <- mgcViz::getViz(mod2)                            # gamViz object of "mod2"

plt_mod_wels <- plot(mgcViz::sm(viz_mod, select = 3))       # select smooth for wels
plt_mod2_wels <- plot(mgcViz::sm(viz_mod2, select = 3))     # select smooth for wels
```

Use the `gridPrint()` function to combine the plots and further fit them as _ggplot_ objects.

``` r
library(ggpubr)
library(grid)
library(gridExtra)

p <- mgcViz::gridPrint(plt_mod_wels +
                           l_fitLine() +
                           l_ciLine() + th2,
                       plt_mod2_wels +
                           l_fitLine() +
                           l_ciLine() + th2, ncol=2)

y_grob <- textGrob("smooth trends",
                   gp=gpar(fontface="plain",
                   col="black", fontsize=18), rot=90)

x_grob <- textGrob("Time",
                   gp=gpar(fontface="plain",
                   col="black", fontsize=18))

top_grob <- textGrob(expression(bold("Smooth effects of wels from two models")),
                                rot = 0, gp = gpar(col = "black", fontsize = 18))

grid.arrange(arrangeGrob(p, left = y_grob, bottom = x_grob, top = top_grob))
```

<p align="center"><img src="/Plots/Fig.sm_eff_5.png?raw=true"></p>

To create a conditional density plot of the partial residuals, we use the parameter _`l_dens(type = "cond")`_.

Let us use the already created object _gamViz_ to plot only the smooth trend of walleye and catfish without the smooth random effects.

``` r
plt_mod_dens <- mgcViz::plot.gamViz(viz_mod, select = c(2,3)) +          # select smooth of pikeperch and wels
                                             l_dens(type = "cond") +
                                             l_fitLine(linetype = 1) +
                                             l_ciLine()
print(plt_mod_dens, pages = 1)
```

<p align="center"><img src="/Plots/Fig.Ex_A1.png?raw=true"></p>


## <a name="heading--4">4. Other packages for GAM visualization: itsadug

[:books:](https://cran.r-project.org/web/packages/itsadug/index.html)`library(itsadug)`

The function `plot_smooth()` is used plot smooth effects on a base R graphs format.

``` r
par(mfrow=c(1,1))
plt_sm_est <- plot_smooth(mod, view = "time",
                               plot_all = c("species"),   # draw fittedlines (e.g. group predictor)
#                              cond=list(species="pike"),
                               sim.ci = TRUE,             # simulated simultaneous CIs (G. Simpson)
                               rm.ranef = FALSE,          # exclude random effects
                               eegAxis = FALSE,           # reverse y-axis
                               hide.label = FALSE,        # fitted values label on the bottom right
#                              plot_all = "fishID",       # draw fitted lines for each individual
#                              legend_plot_all = TRUE,    # disable to show individuals legend
#                              xaxt="n"
                               rug =FALSE,
                               print.summary=FALSE)       # print pair diffs + CIs
# axis(1, col.axis = "transparent", lwd = 1)              # draw x-axis (e.g., if xaxt="n")
```
<p align="center"><img src="/Plots/Fig.sm_eff_6.png?raw=true"></p>

Use the function `plot_smooth()` to plot the random effects smooths.

``` r
par(mfrow=c(1,3), cex=1.1)
plot_smooth(mod, view="time", cond=list(species="pike"), col=cbPalette[1],
                 ylim = c(250, 2500), ylab = "Horizontal range [m]", main='pike')
legend('topleft', legend=c("random-effects.","no random-effects"), # legend settings
                 col=c(cbPalette[1], cbPalette[6]), lwd=2, bty='n')
plot_smooth(mod, view="time", cond=list(species="pike"), col=cbPalette[6],
                 ylim = c(250, 2500), add=TRUE, rm.ranef=TRUE, xpd=TRUE)
plot_smooth(mod, view="time", cond=list(species="pikeperch"), col=cbPalette[1],
                 ylim = c(250, 2500), ylab = "",  main = 'pikeperch')
legend('topleft', legend=c("random-effects.","no random-effects"), # legend settings
                 col=c(cbPalette[1], cbPalette[6]), lwd=2, bty='n')
plot_smooth(mod, view="time", cond=list(species="pikeperch"), col=cbPalette[6],
                 ylim = c(250, 2500), add=TRUE, rm.ranef = TRUE, xpd=TRUE)
plot_smooth(mod, view ="time", cond=list(species="wels"), col=cbPalette[1],
                 ylim = c(250, 2500), ylab = "", main='wels')
plot_smooth(mod, view="time", cond=list(species="wels"), col=cbPalette[6],
                 ylim = c(250, 2500), add=TRUE, rm.ranef=TRUE, xpd=TRUE)
legend('topleft', legend=c("random-effects.","no random-effects"), # legend settings
                 col=c(cbPalette[1], cbPalette[6]), lwd=2, bty='n')
```
<p align="center"><img src="/Plots/Fig.sm_eff_7.png?raw=true"></p>


The functions `plot_diff()` and `plot_diff2()` can be used to show the differences in smoothing trends based on model predictions.

- In this case, we need to make pairwise comparisons and cannot use all groups in the same function. So we end up with three different data sets for comparisons in _ggplot2_

``` r
par(mfrow = c(1,3))
sm_dif_ppk <- plot_diff(mod, view = "time",
                             comp = list(species=c("pike", "pikeperch")),   # pike-pikeperch
                             main = "pike-pikeperch",
                             rm.ranef = TRUE)
sm_dif_pw <-  plot_diff(mod, view="time",
                             comp = list(species = c("pike", "wels")),      # pike-wels
                             main = "pike-wels",
                             rm.ranef = TRUE)
sm_dif_pkw <- plot_diff(mod, view = "time",
                             comp = list(species=c("pikeperch", "wels")),   # pikeperch-wels
                             main = "pikeperch-wels",
                             rm.ranef = TRUE)
```
<p align="center"><img src="/Plots/Fig.sms_diffs_itsa.png?raw=true"></p>


Similarly, for models containing tensor product smooths, we can use the function `itsadug::plot_diff2` to represent differences in smooth surfaces.

The code syntax of this function is similar to `itsadug::plot_diff`, but adds the `view` parameter to specify the two variables that represent the smooth 2D surface.

Also `plot_diff2()` is similar to `mgcv::vis.gam()`, but for 'by' variable pairs the parameter `comp` is used instead of `cond`.

``` r
par(mfrow=c(1,3))
plot_diff2(mod_te, view=c("bl", "temp"),
                   comp=list(species=c("pike", "pikeperch")),
                   main='diff. surface pike-pikeperch',
                   cex.main = 1.8, cex.axis = 1.5,
                   color = "heat", n.grid = 150)
plot_diff2(mod_te, view=c("bl", "temp"),
                   comp=list(species=c("pikeperch", "wels")),
                   main='diff. surface pikeperch-wels',
                   cex.main = 1.8, cex.axis = 1.5,
                   color = "heat", n.grid = 150)
plot_diff2(mod_te, view=c("bl", "temp"),
                   comp=list(species=c("pike", "wels")),
                   main='diff. surface pike-wels',
                   cex.main = 1.8, cex.axis = 1.5,
                   color = "heat", n.grid = 150)
```
```{r, out.width='90%', fig.align='center', fig.cap='Histogram', echo = FALSE}
knitr::include_graphics('Fig.sms_diffs_itsa_tensor.png')
```
<p align="center"><img src="/Plots/Fig.sms_diffs_itsa_tensor.png?raw=true"></p>


## <a name="heading--5">5. Exercises

---
**EXERCISE 1.**

Plot the partial residuals of the model __"mod2"__ , smooth lines in red colour, scaled y-axis with a classic _ggplot_ theme.

<details><summary><strong>Solution</strong></summary><blockquote>

``` r
gratia::draw(mod2, residuals = TRUE,
                   scales = "fixed",
                   smooth_col = "red") &
                   th2
```

<p align="center"><img src="/Plots/Fig.Ex1.png?raw=true"></p>

  </blockquote></details>
</blockquote></details>


---
**EXERCISE 2.**

Plot the smooth time difference between pike and pikeperch from the model __"mod"__ with CIs set to 99%. Use a _ggplot_ theme and mark the significant area of the smooth curve with a different colour.

_Tip 1_: You can filter the species pairs from the tibble table for the smooth difference and adjust the credible intervals of the smooth difference with the parameter `ci_level`.

_Tip 2_: Use `geom_ribbon()` to display upper and lower credible intervals and `stat_smooth` with `method = "gam"` for fitted difference smooth.

<details><summary><strong>Solution</strong></summary><blockquote>

Let us first create a suitable dataframe to work with `ggplot()`.

``` r
sms_diffs_p.pp <- gratia::difference_smooths(mod, smooth = "s(time)",
                                                  ci_level = 0.99)           # credible interval set to 99%

sms_diffs_p.pp$level_1 <- as.factor(sms_diffs_p.pp$level_1)                  # convert the level pair to a factor
sms_diffs_p.pp$level_2 <- as.factor(sms_diffs_p.pp$level_2)

sms_diffs_p.pp <- sms_diffs_p.pp %>%
                    filter(level_1  == "pike" & level_2  == "pikeperch")     # select only the difference smooth 'pike-pikeperch'

sms_diffs_p.pp  <- with(sms_diffs_p.pp,
                        data.frame(diff = as.numeric(diff),
                                  time = as.Date(time),
                                  CI_upper = as.numeric(upper),
                                  CI_lower = as.numeric(lower),
                                  Signif = signifD(sms_diffs_p.pp$diff,
                                                   sms_diffs_p.pp$diff,
                                                   sms_diffs_p.pp$upper,
                                                   sms_diffs_p.pp$lower,
                                                   eval = 0)
                                  )
                        )
```

Finally use `ggplot`to plot the difference smooth displaying the significant period on the fitted curve in a different color.
``` r
ggplot(sms_diffs_p.pp,
       aes(x = time, y = diff)) +
       geom_ribbon(aes(ymax = CI_upper, ymin = CI_lower), alpha = 0.2, fill = cbPalette[4]) +
       stat_smooth(size = 1, method = "gam", level = 0.95, col = cbPalette[4]) +
       geom_line(aes(y = Signif.incr), size = 2,col = cbPalette[3]) +
       geom_line(aes(y = Signif.decr), size = 2,col = cbPalette[2]) +
       scale_x_date(breaks = "1 month", minor_breaks = "1 day", date_labels = "%b %y") +
       geom_hline(data = sms_diffs_p.pp, aes(yintercept = 0), color="black", linetype="dashed") +
       labs(title = "Difference smooth pike-pikeperch", x = "Month", y = "Difference smooth") +
       th2
```

<p align="center"><img src="/Plots/Fig.Ex2.png?raw=true"></p>

  </blockquote></details>
</blockquote></details>



---


***EXERCISE 3***

Draw the derivatives of smooth time for pikeperch using the model __"mod"__ .

_Tip_: You can filter the derivatives for wels from the tibble table and use `geom_ribbon()` to show upper and lower CIs, and `stat_smooth` with _`method = "gam"`_ for the fitted curve.

<details><summary><strong>Solution</strong></summary><blockquote>

Subset the tibble table of derivatives based on pikeperch.

``` r
sms_deriv$smooth <- as.factor(sms_deriv$smooth)                              # 'smooth' to factor
sms_deriv_pp <- sms_deriv %>% filter(smooth  == "s(time):speciespikeperch")  # sel. wels derivatives
sms_deriv_pp  <- with(sms_deriv_pp,
                  data.frame(derivative = as.numeric(derivative),
                             time = as.Date(data),
                             CI_upper = as.numeric(upper),
                             CI_lower = as.numeric(lower)))
```

Plot the derivatives with `ggplot`.
``` r
ggplot(sms_deriv_pp,
       aes(x = time, y = derivative)) +
       geom_hline(data = sms_deriv_pp, aes(yintercept = 0), color="black", linetype="dashed") +
       geom_ribbon(aes(ymax = CI_upper, ymin = CI_lower), alpha = 0.2, fill = cbPalette[4]) +
       stat_smooth(size = 1, method = "gam", level = 0.95, col = cbPalette[4]) +
       scale_x_date(breaks = "1 month", minor_breaks = "1 day", date_labels = "%b %y") +
       ggtitle("Rate of change of wels catfish longitudinal range") +
       ylab(expression(italic(hat(f) * "'") * ("estimated trend"))) + xlab("Month") +
       theme_bw() +
       theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
       theme(plot.title = element_text(size = 16, face = "bold"),
             axis.title.x = element_text(size = 16),
             axis.text.x = element_text(size = 14),
             axis.text.y = element_text(size = 14),
             axis.title.y = element_text(size = 16))
```

<p align="center"><img src="/Plots/Fig.Ex3.png?raw=true"></p>

  </blockquote></details>
</blockquote></details>

---

# <a name="headreferences"></a> References [:page_facing_up:](#headindex)

- **Simpson G (2022)**. _gratia: Graceful ggplot-Based Graphics and Other Functions for GAMs Fitted using mgcv_. R package version 0.7.2., <URL: https://gavinsimpson.github.io/gratia/>.
- **Simpson, G.L.**.  From the bottom of the heap. https://fromthebottomoftheheap.net/.
- **van Rij, J., Wieling, M., Baayen, R., van Rijn, H. (2017)**. itsadug: Interpreting Time Series and Autocorrelated Data Using GAMMs. R package version 2.3
- **Wood, S. N. (2001)**. mgcv: GAMs and Generalized Ridge Regression for R. R News 1(2):20-25
- **Wood, S. N. (2017)**. Generalized additive models: an introduction with R. Second Edition. Boco Raton: CRC Press.


