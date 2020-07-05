fd.est <- function(x, ini.iter = 50, final.iter = 100) {
  mod <- FlexDir::FD.estimation(x, normalize = FALSE, iter.initial.SEM = ini.iter,
                                iter.final.EM = final.iter, verbose = FALSE)
  list(alpha = mod[[ 1 ]], prob = mod[[ 2 ]], tau =  mod[[ 3 ]], loglik =  mod[[ 4 ]])
}
