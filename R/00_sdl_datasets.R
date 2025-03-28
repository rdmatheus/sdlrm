#' @title Stress in Prison
#'
#' @description This data set consists of the stress levels presented by
#'  26 individuals in an experimental study conducted by Verdot et al. (2010)
#'  in a French penitentiary. The inmates were divided into two groups,
#'  one formed by individuals who spontaneously opted to practice sports;
#'  and another one with those who did not wish to perform physical activity.
#'  The observations consist of the stress levels presented by detainees
#'  at the beginning and end of the experiment.
#'
#' @format A data frame with 26 observations on the following 4 variables.
#'
#' \itemize{
#' \item \code{group}: a factor, which identifies whether the individual belongs to the control or the experimental group.
#' \item \code{pss_before}: stress measurement before training.
#' \item \code{pss_after}: stress measurement before training.
#' \item \code{difference}: the difference between the stress levels obtained
#'  at the end of the experiment and at the beginning, that is, \code{pss_after - pss_before}.
#' }
#'
#' @details To measure the stress level, Verdot et al (2010) used the
#' Perceived Stress Scale (PSS) (Cohen, Kamarck and Mermelstein, 1983),
#'  which is a discrete scale and one of the most used psychological
#'  tools to measure the levels of perceived non-specific stress in an
#'  individual.
#'
#' @usage data(pss)
#'
#' @references Cohen, S., Kamarck, T., and Mermelstein, R. (1983). A
#'     global measure of perceived stress. \emph{Journal of Health and Social
#'     Behavior}, \bold{24}, 385---396.
#'
#' @references Verdot, C., Champely, S., Clément, M., and Massarelli, R. (2010).
#'     A simple tool to ameliorate detainees’ mood and well-being in prison: Physical activities.
#'     \emph{International Review on Sport & Violence}, \bold{2}, 83---93.
#'
"pss"



#' @title Relative Renal Function by Scintigraphy
#'
#' @description Data from a study conducted by Domingues et al. (2006) to compare the renal
#'     function measurements of 111 patients with either 99mTc-DTPA or 99mTc-EC dynamic
#'     scintigraphies with that measured using 99mTc-DMSA static scintigraphy. The measurements
#'     reflect the percentage of total renal function in the left kidney and were obtained on a
#'     discrete scale.
#'
#' @format A data frame with 111 observations on the following 4 variables.
#'
#' \itemize{
#' \item \code{static}: measurements obtained by static scintigraphy.
#' \item \code{dynamic}: measurements obtained by dynamic scintigraphy, which may have been
#'      performed with the radiopharmaceutical 99mTc-DTPA or 99mTc-EC.
#' \item \code{difference}: the difference between the renal functions measured using dynamic and
#'     static scintigraphies, that is, \code{dynamic - static}.
#' \item \code{agent}: a factor with the radiopharmaceutical used in the dynamic scintigraphy, with
#'     levels \code{"DTPA"} and \code{"EC"}.
#' \item \code{age}: the patient age, in years.
#' \item \code{sex}: a factor with the patient gender, with levels \code{"F"} and \code{"M"}.
#' }
#'
#' @details Renal scintigraphy is a diagnostic imaging method of nuclear medicine used to measure
#'     kidney function. It is divided between static and dynamic, which differ in the procedure and
#'     the technical evaluation of the results. There are different radiopharmaceuticals used in the
#'     exam. For instance, static renal scintigraphies can use the technetium-99m dimercaptosuccinic
#'     acid (99mTc-DMSA), while dynamic scintigraphies can use the technetium-99m diethylenetriamine
#'     pentaacetic acid (99mTc-DTPA) or the technetium-99m ethylenedicysteine (99mTc-EC). The static
#'     renal agent 99mTc-DMSA is considered the most reliable method for measuring relative renal
#'     function (Kawashima et al., 1998; Martínez et al., 2002). However, this agent has drawbacks, such as
#'     relatively higher radiation dose (Kibar et al., 2003). Its comparison with different
#'     radiopharmaceuticals used in dynamic scintigraphy is of interest in the medical literature.
#'
#' @usage data(scint)
#'
#' @references Domingues, F., Fujikawa, G., Decker, H., Alonso, G., Pereira, J., and Duarte, P. (2006). Comparison
#'    of relative renal function measured with either 99mTc-DTPA or 99mTc-EC dynamic scintigraphies
#'    with that measured with 99mTc-DMSA static scintigraphy. \emph{International Brazilian Journal of
#'    Urology}, \bold{32}, 405---409.
#'
#'    Kibar, M., Yapar, Z., Noyan, A., and Anarat, A. (2003). Technetium-99m-N, N-ethylenedicysteine
#'       and Tc-99m DMSA scintigraphy in the evaluation of renal parenchymal abnormalities in
#'       children. \emph{Annals of Nuclear Medicine}, \bold{17}, 219---225.
#'
#'    Kawashima, A., Sandler, C. M., and Goldman, S. M. (1998). Current roles and controversies in
#'        the imaging evaluation of acute renal infection. \emph{World Journal of Urology},
#'        \bold{16}, 9---17.
#'
#'
#'   Martínez, M., JM, G. D., FJ, D. V., et al. (2002). Comparative study of differential renal
#'       function by DMSA and MAG-3 in congenital unilateral uropathies. \emph{Cirugía Pediátrica},
#'       \bold{15}, 118---121.
#'
#'
"scint"
