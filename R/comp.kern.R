comp.kern <- function(x, type = "alr", h = NULL, thumb = "silverman") {
  if (type == "alr") {  ## additive log-ratio transformation
    y <- Compositional::alr(x)
  } else  y <- alfa(x, 0)$aff
  Compositional::mkde(y, h = h, thumb = thumb)
}
