# Create Channel and subunit classes (V in mV)
setClass("Channel", slots=list(g="numeric", E="numeric"))
setClass("Subunit", slots=list(alpha="function", beta="function"))
setGeneric("tau", function(object, V) standardGeneric("tau"))
setGeneric("infty", function(object, V) standardGeneric("infty"))
setMethod("tau",
          "Subunit",
          function(object, V) {
            return(1/(object@alpha(V)+object@beta(V)))
          })
setMethod("infty",
          "Subunit",
          function(object, V) {
            return(object@alpha(V)/(object@alpha(V)+object@beta(V)))
          })

# Nernst equation
nernst <- function(z, cons_out, cons_in) {
  # Constants
  R <- 8.3144598 # J/mol*K
  F <- 96.485 # J/mol*mV
  T <- 296.15 # Room temperature in Kelvin
  e_chg = 1.602176634e-19 # Elementary charge in Coulombs
  return((R*T)/(z*F)*log(cons_out/cons_in))
}

# Define sub-units and channels
K_orig <- new("Channel", g=36, E=nernst(1,5.4,143))
Na_orig <- new("Channel", g=120, E=nernst(1,145,9.6))
Cl_orig <- new("Channel", g=0.3, E=nernst(-1,127.7,8.7))
Leak_orig <- new("Channel", g=0.3, E=(nernst(-1,127.7,8.7)+nernst(2,1.8,3.5)) )