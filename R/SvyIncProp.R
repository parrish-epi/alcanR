#' Function for calculating the Incidence Proportion
#'
#' @param SurveyDesign Required - a svydesign object created using the survey or srvyr packages.
#' @param event Required - the name of the event outcome, coded as 0 = no event, 1 = event.
#' @param eventime Required - decimal time to end point at time (t). This is decimal age in days from birth to
#'                         event outcome (e.g., first report to OCS), censor (e.g., death or moved out of state),
#'                         or administrative censoring (e.g., end of observation period).
#' @param bystrata1 Optional - strata to calculate the incidence proportion by (e.g., by sex).
#' @param bystrata2 Optional - Second strata classification to calculate the incidence proportion
#'                              by (e.g., by sex and race).
#' @param nreps Default = 100. Number of bootstrap replications. Increase for small values if
#'     needed and decrease for large if taking a long time to run.#'
#' @description This method for calculating the se is based on the approach suggested by
#'     Dr. Lumley using replicate weights. The reason is that the method using ci == TRUE within svykm()
#'     uses an extensive amount of memory for "larger" data sets. This solves that issue and
#'     returns estimates that are similar to his original method. If comparing to the original method it is
#'     important to note that Lumley uses an Aalen (hazard-based) estimator not the Kaplan-Meier estimator.
#'     The implications are that ties and small numbers will result in slight variations in estimates.
#' @details
#' - The standard error of the incidence proportion is calculated using replicate weights on the log(survival) scale.
#' - The method is based on Dr. Lumley's subbootstrap approach, providing estimates similar to his original hazard-based method.
#' - The function handles -Inf or NaN values in replicates by excluding them, but caution is needed if a large proportion of replicates are excluded.
#' - Empirical estimates (mean, median, confidence intervals) from the replicates can become highly variable, especially with small samples.
#' - For subsets (e.g., by groups), the replicate design is derived from the full design object, except when analyzing specific years.
#' @author Jared Parrish, PhD; modified by Richard Trantow (changed svykm.fit to svykm_fit)
#' @import dplyr survey survival
#' @return A data frame containing the incidence proportions, standard errors, confidence intervals, and other summary statistics.
#' @export
svyIncProp <- function(SurveyDesign = NULL, event = NULL,eventime = NULL,bystrata1 = NULL,bystrata2 = NULL,nreps = 100) {
  # require(dplyr)
  # require(survey)
  # require(survival)

  # if ("varlog" %in% names(KM_event)){
  #   stop("oops it looks like you set ci == TRUE in the svykmlist object. Please set ci = FALSE and try again.")
  # } else {

  s1 <- SurveyDesign
  s2 <- s1[order(s1$variables[[eventime]])]

# Overall IP
  if(is.null(bystrata1) & is.null(bystrata2)){
    # Calculate survival object
    SvyKm <- survey::svykm(Surv(get(eventime), get(event))~1, se = F, design = s2)
    SvyDt <- as.data.frame(lapply(SvyKm,identity)) #return list components

    #Update names for scoped variables in design argument
    s2 <- update(s2, died = get(event)==1)
    s2 <- update(s2, time = get(eventime))
    #browser()
    #s1 <- s1[order(s1$variables$time)]


    #Create replicate weights to estimate SE.
    set.seed(3601)
    rpbc <- survey::as.svrepdesign(s2, type="subbootstrap", replicates=nreps, mse=TRUE) #mse=TRUE estimates the SE around the point estimate.
    TIMES <- SvyKm$time                                                         # for estimating survival at each event time
    # Function for calculating the KM object for each replicate
    km_eval_fac<-function(w,data){
      outcome<-with(data, Surv(time, died))
      object<-survey:::svykm_fit(outcome, w)
      idx <- sapply(TIMES, function(t) max(which(object$time <= t), na.rm = TRUE))  # note this results in using the max of the age level for the group with ties for full KM curve
      log(object$surv[idx]) # log survival scale
    }

    sv_vara_reps <- survey::withReplicates(rpbc, km_eval_fac, return.replicates = TRUE)

    # Calculate empirical estimates from replicates
    repInt_rep <- as.data.frame(sv_vara_reps$replicates) %>%
      dplyr::summarise(across(.cols = everything(), ~ 1- exp(quantile(.x, probs = c(0.025,0.500,0.975),na.rm=TRUE)))) %>%
      t() %>%
      as.data.frame() %>%
      select(V2,V3,V1) %>%
      dplyr::rename(empirical.median = V2, empirical.lowerCI = V3, empirical.upperCI = V1)

    repInt_mn_rep <- as.data.frame(sv_vara_reps$replicates) %>%
      #summarise(across(.cols = everything(), ~ 1- quantile(.x, probs = c(0.025,0.500,0.975),na.rm=TRUE))) %>%
      dplyr::summarise(across(.cols = everything(), ~ 1- exp(mean(.x[!is.infinite(.x)], na.rm = TRUE)))) %>%
      t() %>%
      as.data.frame() %>%
      dplyr::rename(empirical.mean = V1)

    # Create bin groups for
    bin.age <- round(max(s2$variables[[eventime]]),0)


    # create Incidence proportion (Ft) at each event time with Confidence intervals

    KM_sub <- cbind(SvyDt,theta = sv_vara_reps$theta, se = SE(sv_vara_reps))

    KM_sub <- KM_sub %>%
      dplyr::group_by(time) %>%
      #filter(surv == max(surv)) %>%
      mutate(time.int = cut(time,breaks=c(0:bin.age),right=FALSE),
             Ft.orig = 1-surv,
             Ft.repl = (1-(exp(theta))),
             Ft.lowerCI = pmax(0,(1-(exp((theta)+(1.96*(se)))))),
             Ft.upperCI = pmin(1,(1-(exp((theta)-(1.96*(se))))))
      ) %>%
      as.data.frame()

    KM_sub <- cbind(KM_sub, repInt_mn_rep, repInt_rep)
    # %>%
    # slice(which.max(exp(time)))

    return(KM_sub)
  }

# Single variable IP by levels of that variable
  if(!is.null(bystrata1) & is.null(bystrata2)){
    # Calculate survival object
    SvyKm <- survey::svykm(Surv(get(eventime), get(event))~strata(get(bystrata1)), se = F, design = s2)
    SvyDt <- do.call(rbind, lapply(names(SvyKm),
                                   function(x){data.frame("strata" = sub(".*=","", x), # only return level, remove variable name
                                                          "time" = (SvyKm[[x]]$time),  # time for every event time
                                                          "surv" = (SvyKm[[x]]$surv)   # surviorship function
                                   )
                                   }
    ))
    #Update names for scoped variables in design argument
    s2 <- update(s2, died = get(event)==1)
    s2 <- update(s2, time = get(eventime))


    # Create replicate object weights to estimate SE.
    # 1 make an empty dataframe
    mydta <- data.frame("strata" = character(),
                        "time" = numeric(),
                        "theta" = numeric(),
                        "se" = numeric(),
                        "empirical.mean" = numeric(),
                        "empirical.median" = numeric(),
                        "empirical.lowerCI" = numeric(),
                        "empirical.upperCI" = numeric())

    # 2 Create replicate design based on survey object
    set.seed(3601)
    TempD <- survey::as.svrepdesign(s2, type="subbootstrap", replicates=nreps, mse=TRUE)

    # 3 create bin ages for estimates
    # Create bin groups for
    bin.age <- round(max(s2$variables[[eventime]]),0)
    # 4 Loop across each strata, create replicates and bootstrap together
    # Note: With the sapply, estimation of S(t) occurs at the max of the object time
    #       this results in the survival/IP curve starting at the first event time opposed to 0/1.

    for (j in levels(as.factor((s2$variables[[bystrata1]])))) {   # Loop across each strata level

      # Subset replicate design for each strata subset
      subTempD <- subset(TempD, get(bystrata1) == j)   # NOTE: This is the best method opposed to creating a rep design on subset design object
      ## note will need to un-comment these parts if adding back in the time points specification.
      #if (is.null(time point)) {
      # Create time variable for age at each event time
      TIMES1 <- SvyDt %>% dplyr::filter(strata %in% j) %>% select(time)
      TIMES <- TIMES1$time
      #} else {
      ###Note: if a time point is specified this is the IP at the point in time and does not include the interval after.
      #     TIMES <- time point
      #   }

      # Function for calculating the KM object for each replicate
      km_eval_fac<-function(w,data){
        outcome<-with(data, Surv(time, died))
        object<-survey:::svykm_fit(outcome, w)
        idx <- sapply(TIMES, function(t) max(which(object$time <= t)))  # note this results in using the max of the age level for the group with ties for full KM curve
        log(object$surv[idx])
      }
      # bootstrap replicates together as a data frame
      #sv_var_reps <- as.data.frame(withReplicates(subTempD, km_eval_fac))
      sv_vara_reps <- survey::withReplicates(subTempD, km_eval_fac, return.replicates = TRUE)

      # Calculate empirical estimates from replicates
      repInt_rep <- as.data.frame(sv_vara_reps$replicates) %>%
        dplyr::summarise(across(.cols = everything(), ~ 1- exp(quantile(.x, probs = c(0.025,0.500,0.975),na.rm=TRUE)))) %>%
        t() %>%
        as.data.frame() %>%
        dplyr::select(V2,V3,V1) %>%
        dplyr::rename(empirical.median = V2, empirical.lowerCI = V3, empirical.upperCI = V1)

      repInt_mn_rep <- as.data.frame(sv_vara_reps$replicates) %>%
        #summarise(across(.cols = everything(), ~ 1- quantile(.x, probs = c(0.025,0.500,0.975),na.rm=TRUE))) %>%
        dplyr::summarise(across(.cols = everything(), ~ 1- exp(mean(.x[!is.infinite(.x)], na.rm = TRUE)))) %>%
        t() %>%
        as.data.frame() %>%
        dplyr::rename(empirical.mean = V1)

      # create Incidence proportion (Ft) at each event time with Confidence intervals
      SvyDt_sub <- SvyDt %>% dplyr::filter(strata %in% j)

      tmpdta <- cbind(SvyDt_sub, theta = sv_vara_reps$theta,
                      se = SE(sv_vara_reps), repInt_mn_rep, repInt_rep)

      # sv_var$strata <- j  # create column specifying strata level
      # sv_var$time <- TIMES # create column specifying the time levels
      mydta <- rbind(mydta,tmpdta) # bind all strata together into one data frame
    }
    # 5 create Ft by group
    KM_sub <- mydta %>% dplyr::group_by(strata,time) %>%
      #slice(which.max(exp(time))) %>%
      mutate(time.int = cut(time,breaks=c(0:bin.age),right=FALSE),
             Ft.orig = 1-surv,
             Ft.repl = (1-(exp(theta))),
             Ft.lowerCI = pmax(0,(1-(exp((theta)+(1.96*(se)))))),
             Ft.upperCI = pmin(1,(1-(exp((theta)-(1.96*(se))))))
      ) %>%
      as.data.frame() %>%
      dplyr::rename(!!sym(bystrata1) := strata)
    return(KM_sub)
  }

# Bivariate  IP by levels of the variable combinations
  if(!is.null(bystrata1) & !is.null(bystrata2)){
    # Calculate survival object
    SvyKm <- survey::svykm(Surv(get(eventime), get(event))~strata(get(bystrata1))+strata(get(bystrata2)), se = F, design = s2)

    st1 <- is.numeric(s2$variables[[bystrata1]])

    SvyDt <- do.call(rbind, lapply(names(SvyKm),
                                   function(x){data.frame("strata1" = if(st1 == TRUE){
                                     gsub(".*[=]([^.]+)[.].*","\\1",x) # only return level, remove variable name}
                                   } else {
                                     sub("[.].*","",x)
                                   },
                                   "strata2" = sub(".*[.|=]","", x), # return level after var name
                                   "time" = (SvyKm[[x]]$time),  # time for every event time
                                   "surv" = (SvyKm[[x]]$surv)   # surviorship function
                                   )
                                   }
    ))
    #Update names for scoped variables in design argument
    s2 <- update(s2, died = get(event)==1)
    s2 <- update(s2, time = get(eventime))


    #Create replicate object weights to estimate SE.

    # 1 make an empty dataframe
    mydta <- data.frame("strata1" = character(),
                        "strata2" = character(),
                        "time" = numeric(),
                        "theta" = numeric(),
                        "se" = numeric(),
                        "empirical.mean" = numeric(),
                        "empirical.median" = numeric(),
                        "empirical.lowerCI" = numeric(),
                        "empirical.upperCI" = numeric())

    # 2 Create replicate design based on survey object
    set.seed(3601)
    TempD <- survey::as.svrepdesign(s2, type="subbootstrap", replicates=nreps, mse=TRUE)
    # 3 create bin ages for estimates
    # Create bin groups for
    bin.age <- round(max(s2$variables[[eventime]]),0)

    # 4 Loop across each strata, create replicates and bootstrap together
    # Note: One minor issue: with the sapply, estimation of S(t) occurs at the max of the object time
    #       this results in the survival/IP curve starting at the first event time opposed to 0/1.

    for (j in levels(as.factor(s2$variables[[bystrata1]]))) {     # Loop across each strata1 level
      for (k in levels(as.factor(s2$variables[[bystrata2]]))){      # Loop across each strata2 level
        # Subsect replicate design by each strata subset

        SubTempD <- subset(TempD, get(bystrata1) == j & get(bystrata2) == k)
        ## note will need to un-comment these parts if adding back in the timepoints specification.
        #if (is.null(timepoint)) {
        # Create time variable for age at each event time
        TIMES1 <- SvyDt %>% dplyr::filter(strata1 %in% j & strata2 %in% k) %>% select(time)
        TIMES <- TIMES1$time
        #} else {
        ###Note: if a timepoint is specified this is the IP at the point in time and does not include the interval after.
        #     TIMES <- timepoint
        #   }

        # Function for calcuating the KM object for each replicate
        km_eval_fac<-function(w,data){
          outcome<-with(data, Surv(time, died))
          object<-survey:::svykm_fit(outcome, w)
          idx <- sapply(TIMES, function(t) max(which(object$time <= t)))  # note this results in using the max of the age level for the group with ties for full KM curve
          log(object$surv[idx])
        }
        # bootstrap replicates together as a dataframe

        sv_vara_reps <- survey::withReplicates(SubTempD, km_eval_fac, return.replicates = TRUE)

        #m[!is.infinite(m)]
        # Calculate empirical estimates from replicates
        repInt_rep <- as.data.frame(sv_vara_reps$replicates) %>%
          dplyr::summarise(across(.cols = everything(), ~ 1- exp(quantile(.x, probs = c(0.025,0.500,0.975),na.rm=TRUE)))) %>%
          t() %>%
          as.data.frame() %>%
          dplyr::select(V2,V3,V1) %>%
          dplyr::rename(empirical.median = V2, empirical.lowerCI = V3, empirical.upperCI = V1)

        repInt_mn_rep <- as.data.frame(sv_vara_reps$replicates) %>%
          #summarise(across(.cols = everything(), ~ 1- quantile(.x, probs = c(0.025,0.500,0.975),na.rm=TRUE))) %>%
          dplyr::summarise(across(.cols = everything(), ~ 1- exp(mean(.x[!is.infinite(.x)], na.rm = TRUE)))) %>%
          t() %>%
          as.data.frame() %>%
          dplyr::rename(empirical.mean = V1)
        #browser()
        # create data set for each strata at each event time with Confidence intervals
        SvyDt_sub <- SvyDt %>% dplyr::filter(strata1 %in% j & strata2 %in% k)
        tmpdta <- cbind(SvyDt_sub, theta = sv_vara_reps$theta,
                        se = SE(sv_vara_reps), repInt_mn_rep, repInt_rep)

        # sv_var$strata <- j  # create column specifying strata level
        # sv_var$time <- TIMES # create column specifying the time levels
        mydta <- rbind(mydta,tmpdta) # bind all strata together into one data frame

        #mydta <- rbind(mydta,sv_var) # bind all strata together into one dataframe
      }
    }

    KM_sub <- mydta %>% dplyr::group_by(strata1,strata2,time) %>%
      #slice(which.max(exp(time))) %>%
      mutate(time.int = cut(time,breaks=c(0:bin.age),right=FALSE),
             Ft.orig = 1-surv,
             Ft.repl = (1-(exp(theta))),
             Ft.lowerCI = pmax(0,(1-(exp((theta)+(1.96*(se)))))),
             Ft.upperCI = pmin(1,(1-(exp((theta)-(1.96*(se))))))
      ) %>%
      as.data.frame() %>%
      dplyr::rename(!!sym(bystrata1) := strata1,
             !!sym(bystrata2) := strata2)
    return(KM_sub)
  }

}

