library(tercen)
library(dplyr)
library(ggplot2)
library(drc)

drcFit <- function(df, x_multiplier = 1) {
  dataY     <- df$.y
  drc_fit   <- try(drm(.y ~ .x, data = df, fct = LL.4(), logDose = 10, na.action=na.omit), silent=TRUE)
  fit_error <- FALSE
  xfit      <- df$.x
  if (class(drc_fit)[1] == 'try-error') {
    fit_error <- TRUE
  } else {
    cfs         <- as.vector(coef(drc_fit))
    hillslope   <- cfs[1]
    Ymin        <- cfs[2]
    Ymax        <- cfs[3]
    MI          <- log10((Ymin/Ymax)^(sign(hillslope))) 
    EC50        <- cfs[4]
    logEC50     <- log10(EC50)
    Rsq         <- 1-(sum(residuals(drc_fit)^2)/sum((dataY-mean(dataY))^2))
    lpRes       <- -log10(var.test(residuals(drc_fit), dataY - mean(dataY))$p.value)
    yfit        <- predict(drc_fit)
    
    err.95  <- try(as.vector((coef(drc_fit) - confint(drc_fit)) [,1]), silent = TRUE)
    if (class(err.95)[1] == 'try-error') {
      fit_error <- TRUE
    }
    else {
      err.95.hillslope = err.95[1]
      err.95.Ymin      = err.95[2]
      err_95.Ymax      = err.95[3]
      err.95.EC50      = err.95[4]
    }
  }
  if (fit_error) {
    hillslope <- Ymin <- Ymax <- EC50 <- logEC50 <- NaN
    MI        <- Rsq  <- lpRes <- 0
    err.95.hillslope <- err.95.Ymin <- err_95.Ymax <- err.95.EC50 <- NaN
    yfit <- df$.y
  }
  
  data.frame(.ri              = df$.ri[1], 
             .ci              = df$.ci[1], 
             Ymax             = Ymax, 
             err_95.Ymax      = err_95.Ymax,
             Ymin             = Ymin, 
             err.95.Ymin      = err.95.Ymin,           
             EC50             = x_multiplier * EC50, 
             err.95.EC50      = x_multiplier * err.95.EC50,
             logEC50          = logEC50, 
             hillslope        = hillslope, 
             err.95.hillslope = err.95.hillslope,
             R2               = Rsq,
             lpRes            = lpRes,
             MI               = MI,
             xfit             = xfit,
             yfit             = yfit
  )
}

ctx = tercenCtx()

x_axis_multiplier <- ifelse(is.null(ctx$op.value('X-axis multiplier')), 1.0, as.double(ctx$op.value('X-axis multiplier')))

ctx %>% 
  dplyr::select(.ri, .ci, .y, .x) %>%
  arrange(.ri, .ci, .x) %>%
  group_by(.ri, .ci) %>%
  do(drcFit(., x_axis_multiplier)) %>%
  ungroup() %>%
  mutate_at(vars(setdiff(colnames(data), c(".ri", ".ci"))), ~ifelse(!is.nan(.) & abs(.) > 1e30, NaN, .)) %>%
  ctx$addNamespace() %>%
  ctx$save()
