#' @title Stress in prison
#'
#' @description This data set consists of the stress levels presented by
#'  26 individuals in an experimental study conducted by Verdot et al (2008)
#'  in a French penitentiary. The inmates were divided into two groups,
#'  one formed by individuals who spontaneously opted to practice sports;
#'  and another one with those who did not wish to perform physical activity.
#'  The observations consist of the stress levels presented by detainees
#'  at the beginning and end of the experiment.
#'
#' @format A data frame with 26 observations on the following 4 variables.
#'
#' \itemize{
#' \item \code{Group}: a factor, which identifies whether the individual belongs to the control or the experimental group.
#' \item \code{PSSbefore}: stress measurement before training.
#' \item \code{PSSafter}: stress measurement before training.
#' \item \code{Diff}: the difference between the stress levels obtained
#'  at the end of the experiment and at the beginning: \code{PSSafter} - \code{PSSbefore}.
#' }
#'
#' @details To measure the stress level, Verdot et al (2008) used the
#' Perceived Stress Scale (PSS) (Cohen, Kamarck and Mermelstein, 1983),
#'  which is a discrete scale and one of the most used psychological
#'  tools to measure the levels of perceived non-specific stress in an
#'  individual.
#'
#' @usage data(pss)
#'
#' @references Cohen, S., Kamarck, T., and Mermelstein, R. (1983). A
#'     global measure of perceived stress. Journal of Health and Social
#'     Nehavior, 24, 385-396.
#'
#' @references Verdot, C., Champely, S., Massarelli, R., and Clement,
#'     M. (2008). Physical activities in prison as a tool to ameliorate
#'     detainees mood and well-being. International Review on Sport
#'     Violence: Sport in prison, 2, 83-93.
#'
#' @references Medeiros, R. M. R., Bourguignon, M. (2020). Mean and
#'     dispersion regression model for paired count data with
#'     application to quality of life in prisons.
#'
"pss"
