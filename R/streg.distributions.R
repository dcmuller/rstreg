streg.distributions <- vector(mode="list")
streg.distributions[["weibull"]] <- list(
  name="Weibull",
  lh = function(xb, lp, tt) {xb + exp(lp)*log(tt)},
  dlh = function(x,z,lp,tt) {
    cbind(x,
          z*exp(lp)*log(tt))
  },
  lhr = function(xb, lp, xb2, lp2, tt) {xb + exp(lp)*log(tt) - xb2 - exp(lp2)*log(tt)},
  dlhr = function(x, z, x2, z2, lp, lp2, tt) {
    cbind(x-x2,
          z * exp(lp)*log(tt) - z2 * exp(lp2)*log(tt))
  },
  lH = function(xb,lp,tt,tt0) {
    xb + log(exp(exp(lp)*log(tt)) - ifelse(tt0>0, exp(exp(lp)*log(tt0)), 0))
  },
  dlH = function(x,z,lp,tt,tt0) {
    cbind(x,
          z * ((exp(lp) *
                  (log(tt)*exp(exp(lp)*log(tt)) - ifelse(tt0>0, log(tt0)*exp(exp(lp)*log(tt0)),0)))/
                 ((exp(exp(lp)*log(tt)) - ifelse(tt0>0, exp(exp(lp)*log(tt0)),0)))) )
  },
  qf = function(qq,xb,lp) {exp(1/exp(lp)*log(-log((1-qq))/exp(xb)))},
  dqf = function(qq,x,z,xb,lp) {
    cbind(
      -(x * (exp(log(-(log(1 - qq)/exp(xb)))/exp(lp))/exp(lp))),
      -(z * (exp(log(-(log(1 - qq)/exp(xb)))/exp(lp)) * (log(-(log(1 - qq)/exp(xb)))/exp(lp))))
    )
  }
)
