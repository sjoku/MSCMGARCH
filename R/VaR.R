# Function to calculate the Value at Risk (VaR) for a Gaussian distribution
VaR_Gaussian <- function(portfolio_weights, H, level) {
  SD = portfolio_weights %*% H %*% portfolio_weights
  VaR = qnorm(level, mean = 0, sd = sqrt(SD))
  return(VaR)
}

# Function to calculate the Expected Shortfall (ES) for a Gaussian distribution
ES_Gaussian <- function(portfolio_weights, H, level, VaR) {
  SD = portfolio_weights %*% H %*% portfolio_weights
  ES = (1 / level) * integrate(function(z) { z * dnorm(z, mean = 0, sd = sqrt(SD)) }, lower = -Inf, upper = VaR, rel.tol = 1e-15)$value
  return(ES)
}

# Function to calculate the u statistic for a Gaussian distribution
u_Gaussian <- function(portfolio_weights, H, level, VaR, portfolio_return) {
  SD = portfolio_weights %*% H %*% portfolio_weights
  u = 0
  if (portfolio_return <= VaR) {
    u = pnorm(portfolio_return, mean = 0, sd = sqrt(SD))
  }
  return(u)
}

# Function to calculate the Value at Risk (VaR) for a CMGARCH model
VaR_CMGARCH <- function(portfolio_weights, H, level, type, par) {
  correlation_matrix = cor_mat(par, type)
  portfolio = portfolio_weights %*% eigen_value_decomposition(H) %*% solve(t(chol(correlation_matrix)))
  if (portfolio[2] == 0 && portfolio[3] == 0) {
    if (portfolio[1] > 0) {
      VaR = qnorm(level, sd = portfolio[1])
    } else {
      VaR = 1 - qnorm(level, sd = portfolio[1])
    }
    return(VaR)
  }
  if (portfolio[1] == 0 && portfolio[2] == 0) {
    if (portfolio[3] > 0) {
      VaR = qnorm(level, sd = portfolio[3])
    } else {
      VaR = 1 - qnorm(level, sd = portfolio[3])
    }
    return(VaR)
  }
  if (portfolio[1] == 0 && portfolio[3] == 0) {
    if (portfolio[2] > 0) {
      VaR = qnorm(level, sd = portfolio[2])
    } else {
      VaR = 1 - qnorm(level, sd = portfolio[2])
    }
    return(VaR)
  }
  if (portfolio[3] == 0) {
    f = function(z) {
      if (portfolio[1] >= 0) {
        return(cubintegrate(function(u) { copulaH2(pnorm(z / portfolio[1] - portfolio[2] * u / portfolio[1]), pnorm(u), par[1], type[1]) * dnorm(u) }, lower = -2, upper = 2)$integral - level)
      } else {
        return((1 - cubintegrate(function(u) { copulaH2(pnorm(z / portfolio[1] - portfolio[2] * u / portfolio[1]), pnorm(u), par[1], type[1]) * dnorm(u) }, lower = -2, upper = 2)$integral) - level)
      }
    }
    VaR = uniroot(f, interval = c(-1, 0))$root
    return(VaR)
  }
  if (portfolio[2] == 0) {
    f = function(z) {
      if (portfolio[1] >= 0) {
        return(cubintegrate(function(u) { copulaH2(pnorm(z / portfolio[1] - portfolio[3] / portfolio[1]), pnorm(u), par[2], type[2]) * dnorm(u) }, lower = -2, upper = 2)$integral - level)
      } else {
        return((1 - cubintegrate(function(u) { copulaH2(pnorm(z / portfolio[1] - portfolio[3] / portfolio[1]), pnorm(u), par[2], type[2]) * dnorm(u) }, lower = -2, upper = 2)$integral) - level)
      }
    }
    VaR = uniroot(f, interval = c(-1, 0))$root
    return(VaR)
  }
  if (portfolio[1] == 0) {
    f = function(z) {
      if (portfolio[2] >= 0) {
        return(cubintegrate(function(u) { copulaH2(pnorm(z / portfolio[2] - u * portfolio[3] / portfolio[2]), pnorm(u), par[3], type[3]) * dnorm(u) }, lower = -2, upper = 2)$integral - level)
      } else {
        return((1 - cubintegrate(function(u) { copulaH2(pnorm(z / portfolio[2] - u * portfolio[3] / portfolio[2]), pnorm(u), par[3], type[3]) * dnorm(u) }, lower = -2, upper = 2)$integral) - level)
      }
    }
    VaR = uniroot(f, interval = c(-1, 0))$root
    return(VaR)
  } else {
    if (portfolio[1] > 0) {
      f = function(k) { integral(function(x, y) { VineH2(pnorm(k / portfolio[1] - x * portfolio[2] / portfolio[1] - y * portfolio[3] / portfolio[1]), pnorm(x), pnorm(y), par, type) * dnorm(x) * dnorm(y) }, bounds = list(x = c(-4, 4), y = c(-4, 4)))$value }
      norm_const = f(-1)
      g = function(k) { f(k) - norm_const - level }
    } else {
      f = function(k) { 1 - integral(function(x, y) { VineH2(pnorm(k / portfolio[1] - x * portfolio[2] / portfolio[1] - y * portfolio[3] / portfolio[1]), pnorm(x), pnorm(y), par, type) * dnorm(x) * dnorm(y) }, bounds = list(x = c(-4, 4), y = c(-4, 4)))$value }
      norm_const = f(-1)
      g = function(k) { f(k) - norm_const - level }
    }
    VaR = uniroot(g, interval = c(-1, 0))$root
    return(VaR)
  }
}

# Function to calculate the Expected Shortfall (ES) for a CMGARCH model
ES_CMGARCH <- function(portfolio_weights, H, level, type, par, VaR) {
  correlation_matrix = cor_mat(par, type)
  portfolio = portfolio_weights %*% eigen_value_decomposition(H) %*% solve(t(chol(correlation_matrix)))
  if (portfolio[2] == 0 && portfolio[3] == 0) {
    ES = (integrate(function(z) { z * dnorm(z, sd = portfolio[1]) }, lower = -Inf, upper = VaR, rel.tol = 1e-15)$value) / level
    return(ES)
  }
  if (portfolio[1] == 0 && portfolio[2] == 0) {
    ES = (integrate(function(z) { z * dnorm(z, sd = portfolio[3]) }, lower = -Inf, upper = VaR, rel.tol = 1e-15)$value) / level
    return(ES)
  }
  if (portfolio[1] == 0 && portfolio[3] == 0) {
    ES = (integrate(function(z) { z * dnorm(z, sd = portfolio[2]) }, lower = -Inf, upper = VaR)$value) / level
    return(ES)
  }
  if (portfolio[3] == 0) {
    f = function(z) {
      if (portfolio[1] > 0) {
        return(cubintegrate(function(u) { copulaH2(pnorm(z / portfolio[1] - u * portfolio[2] / portfolio[1]), pnorm(u), par[1], type[1]) * dnorm(u) }, lower = -2, upper = 2)$integral)
      } else {
        return((1 - cubintegrate(function(u) { copulaH2(pnorm(z / portfolio[1] - u * portfolio[2] / portfolio[1]), pnorm(u), par[1], type[1]) * dnorm(u) }, lower = -2, upper = 2)$integral))
      }
    }
    ES = -integral(f, bounds = list(z = c(-1, VaR)))$value
    return(ES / f(VaR))
  }
  if (portfolio[2] == 0) {
    f = function(z) {
      if (portfolio[1] > 0) {
        return(integrate(function(u) { return(copulaH2(pnorm(z / portfolio[1] - u * portfolio[3] / portfolio[1]), pnorm(u), par[2], type[2]) * dnorm(u)) }, lower = -2, upper = 2)$value)
      } else {
        return((1 - integrate(function(u) { return(copulaH2(pnorm(z / portfolio[1] - u * portfolio[3] / portfolio[1]), pnorm(u), par[2], type[2]) * dnorm(u)) }, lower = -2, upper = 2)$value))
      }
    }
    ES = -integral(f, bounds = list(z = c(-1, VaR)))$value
    return(ES / f(VaR))
  }
  if (portfolio[1] == 0) {
    f = function(z) {
      if (portfolio[2] > 0) {
        return(cubintegrate(function(u) { copulaH2(pnorm(z / portfolio[2] - u * portfolio[3] / portfolio[2]), pnorm(u), par[3], type[3]) * dnorm(u) }, lower = -2, upper = 2)$integral)
      } else {
        return((1 - cubintegrate(function(u) { copulaH2(pnorm(z / portfolio[2] - u * portfolio[3] / portfolio[2]), pnorm(u), par[3], type[3]) * dnorm(u) }, lower = -2, upper = 2)$integral))
      }
    }
    ES = -integral(f, bounds = list(z = c(-1, VaR)))$value
    return(ES / f(VaR))
  } else {
    f = function(x) { integral(function(y, z) { VineCopula(c(pnorm(x / portfolio[1] - y * portfolio[2] / portfolio[1] - z * portfolio[3] / portfolio[1]), pnorm(y), pnorm(z)), as.matrix(par), type) * dnorm(x / portfolio[1] - y * portfolio[2] / portfolio[1] - z * portfolio[3] / portfolio[1]) * dnorm(y) * dnorm(z) }, bounds = list(y = c(-4, 4), z = c(-4, 4)), absTol = 1e-3)$value }
    norm_const = integral(function(z) { f(z) }, bounds = list(z = c(-1, VaR)), absTol = 1e-3)$value
    ES = integral(function(z) { z * (f(z)) }, bounds = list(z = c(-1, VaR)), absTol = 1e-3)$value / norm_const
    return(ES)
  }
}

# Function to calculate the u statistic for a CMGARCH model
u_CMGARCH <- function(portfolio_weights, H, level, type, par, VaR, portfolio_return) {
  u = 0
  if (portfolio_return <= VaR) {
    correlation_matrix = cor_mat(par, type)
    portfolio = portfolio_weights %*% eigen_value_decomposition(H) %*% solve(t(chol(correlation_matrix)))
    if (portfolio[2] == 0 && portfolio[3] == 0) {
      u = pnorm(portfolio_return, sd = portfolio[1])
    }
    if (portfolio[1] == 0 && portfolio[2] == 0) {
      u = pnorm(portfolio_return, sd = portfolio[3])
    }
    if (portfolio[1] == 0 && portfolio[3] == 0) {
      u = pnorm(portfolio_return, sd = portfolio[2])
    }
    if (portfolio[3] == 0) {
      f = function(z) {
        if (portfolio[1] > 0) {
          return(cubintegrate(function(u) { copulaH2(pnorm(z / portfolio[1] - portfolio[2] * u / portfolio[1]), pnorm(u), par[1], type[1]) * dnorm(u) }, lower = -2, upper = 2)$integral)
        } else {
          return(1 - cubintegrate(function(u) { copulaH2(pnorm(z / portfolio[1] - portfolio[2] * u / portfolio[1]), pnorm(u), par[1], type[1]) * dnorm(u) }, lower = -2, upper = 2)$integral)
        }
      }
      u = f(portfolio_return)
    }
    if (portfolio[2] == 0) {
      f = function(z) {
        if (portfolio[1] > 0) {
          return(cubintegrate(function(u) { copulaH2(pnorm(z / portfolio[1] - portfolio[3] * u / portfolio[1]), pnorm(u), par[2], type[2]) * dnorm(u) }, lower = -2, upper = 2)$integral)
        } else {
          return(1 - cubintegrate(function(u) { copulaH2(pnorm(z / portfolio[1] - portfolio[3] * u / portfolio[1]), pnorm(u), par[2], type[2]) * dnorm(u) }, lower = -2, upper = 2)$integral)
        }
      }
      u = f(portfolio_return)
    }
    if (portfolio[1] == 0) {
      f = function(z) {
        if (portfolio[2] > 0) {
          return(cubintegrate(function(u) { copulaH2(pnorm(z / portfolio[2] - u * portfolio[3] / portfolio[2]), pnorm(u), par[3], type[3]) * dnorm(u) }, lower = -2, upper = 2)$integral)
        } else {
          return(1 - cubintegrate(function(u) { copulaH2(pnorm(z / portfolio[2] - u * portfolio[3] / portfolio[2]), pnorm(u), par[3], type[3]) * dnorm(u) }, lower = -2, upper = 2)$integral)
        }
      }
      u = f(portfolio_return)
    } else {
      if (portfolio[1] > 0) {
        f = function(k) { integral(function(x, y) { VineH2(pnorm(k / portfolio[1] - x * portfolio[2] / portfolio[1] - y * portfolio[3] / portfolio[1]), pnorm(x), pnorm(y), par, type) * dnorm(x) * dnorm(y) }, bounds = list(x = c(-4, 4), y = c(-4, 4)))$value }
        norm_const = f(-1)
        g = function(k) { f(k) - norm_const }
      } else {
        f = function(k) { 1 - integral(function(x, y) { VineH2(pnorm(k / portfolio[1] - x * portfolio[2] / portfolio[1] - y * portfolio[3] / portfolio[1]), pnorm(x), pnorm(y), par, type) * dnorm(x) * dnorm(y) }, bounds = list(x = c(-4, 4), y = c(-4, 4)))$value }
        norm_const = f(-1)
        g = function(k) { f(k) - norm_const }
      }
      u = g(portfolio_return)
    }
  }
  return(u)
}

# Function to calculate the Value at Risk (VaR) for a two-regime MSCMGARCH model
VaR_MSCMGARCH <- function(portfolio_weights, H, level, type, par, filterprobs) {
  if (length(type) == 3 && length(par) == 5 && !is.null(filterprobs)) {
    correlation_matrix = cor_mat(par[3:5], type)
    p = par[1]
    q = par[2]
    par = par[3:5]
    SD = sqrt(portfolio_weights %*% H %*% portfolio_weights)
    portfolio = portfolio_weights %*% eigen_value_decomposition(H) %*% solve(t(chol(correlation_matrix)))
    if (sum(portfolio == 0) > 1) {
      g = function(k) {
        if (sum(portfolio) > 0) {
          VaR = (p * filterprobs + (1 - q) * (1 - filterprobs)) * pnorm(k, sd = SD) + ((1 - p) * filterprobs + (q) * (1 - filterprobs)) * pnorm(k, sd = sum(portfolio)) - level
        } else {
          VaR = (p * filterprobs + (1 - q) * (1 - filterprobs)) * pnorm(k, sd = SD) + ((1 - p) * filterprobs + (q) * (1 - filterprobs)) * (1 - pnorm(k, sd = abs(sum(portfolio)))) - level
        }
        return(VaR)
      }
      return(uniroot(g, interval = c(-1, 0))$root)
    }
    if (portfolio[3] == 0) {
      f = function(z) {
        if (portfolio[1] > 0) {
          return((p * filterprobs + (1 - q) * (1 - filterprobs)) * pnorm(z, mean = 0, sd = SD) + ((1 - p) * filterprobs + (q) * (1 - filterprobs)) * cubintegrate(function(u) { copulaH2(pnorm(z / portfolio[1] - u * portfolio[2] / portfolio[1]), pnorm(u), par[1], type[1]) * dnorm(u) }, lower = -2, upper = 2)$integral - level)
        } else {
          return((p * filterprobs + (1 - q) * (1 - filterprobs)) * pnorm(z, mean = 0, sd = SD) + ((1 - p) * filterprobs + (q) * (1 - filterprobs)) * (1 - cubintegrate(function(u) { copulaH2(pnorm(z / portfolio[1] - u * portfolio[2] / portfolio[1]), pnorm(u), par[1], type[1]) * dnorm(u) }, lower = -2, upper = 2)$integral) - level)
        }
      }
      VaR = uniroot(f, interval = c(-1, 0))$root
      return(VaR)
    }
    if (portfolio[2] == 0) {
      f = function(z) {
        if (portfolio[1] > 0) {
          return((p * filterprobs + (1 - q) * (1 - filterprobs)) * pnorm(z, mean = 0, sd = SD) + ((1 - p) * filterprobs + (q) * (1 - filterprobs)) * cubintegrate(function(u) { return(copulaH2(pnorm(z / portfolio[1] - u * portfolio[3] / portfolio[1]), pnorm(u), par[2], type[2]) * dnorm(u)) }, lower = -2, upper = 2)$integral - level)
        } else {
          return((p * filterprobs + (1 - q) * (1 - filterprobs)) * pnorm(z, mean = 0, sd = SD) + ((1 - p) * filterprobs + (q) * (1 - filterprobs)) * (1 - cubintegrate(function(u) { return(copulaH2(pnorm(z / portfolio[1] - u * portfolio[3] / portfolio[1]), pnorm(u), par[2], type[2]) * dnorm(u)) }, lower = -2, upper = 2)$integral) - level)
        }
      }
      VaR = uniroot(f, interval = c(-1, 0))$root
      return(VaR)
    }
    if (portfolio[1] == 0) {
      f = function(z) {
        if (portfolio[2] > 0) {
          return((p * filterprobs + (1 - q) * (1 - filterprobs)) * pnorm(z, sd = SD) + ((1 - p) * filterprobs + (q) * (1 - filterprobs)) * cubintegrate(function(u) { copulaH2(pnorm(z / portfolio[2] - u * portfolio[3] / portfolio[2]), pnorm(u), par[3], type[3]) * dnorm(u) }, lower = -2, upper = 2)$integral - level)
        } else {
          return((p * filterprobs + (1 - q) * (1 - filterprobs)) * pnorm(z, sd = SD) + ((1 - p) * filterprobs + (q) * (1 - filterprobs)) * (1 - cubintegrate(function(u) { copulaH2(pnorm(z / portfolio[2] - u * portfolio[3] / portfolio[2]), pnorm(u), par[3], type[3]) * dnorm(u) }, lower = -2, upper = 2)$integral) - level)
        }
      }
      VaR = uniroot(f, interval = c(-1, 0))$root
      return(VaR)
    } else {
      if (portfolio[1] > 0) {
        f = function(k) { integral(function(x, y) { VineH2(pnorm(k / portfolio[1] - x * portfolio[2] / portfolio[1] - y * portfolio[3] / portfolio[1]), pnorm(x), pnorm(y), par, type) * dnorm(x) * dnorm(y) }, bounds = list(x = c(-4, 4), y = c(-4, 4)))$value }
        norm_const = f(-1)
        g = function(k) { (p * filterprobs + (1 - q) * (1 - filterprobs)) * pnorm(k, sd = SD) + ((1 - p) * filterprobs + (q) * (1 - filterprobs)) * (f(k) - norm_const) - level }
      } else {
        f = function(k) { 1 - integral(function(x, y) { VineH2(pnorm(k / portfolio[1] - x * portfolio[2] / portfolio[1] - y * portfolio[3] / portfolio[1]), pnorm(x), pnorm(y), par, type) * dnorm(x) * dnorm(y) }, bounds = list(x = c(-4, 4), y = c(-4, 4)))$value }
        norm_const = f(-1)
        g = function(k) { (p * filterprobs + (1 - q) * (1 - filterprobs)) * pnorm(k, sd = SD) + ((1 - p) * filterprobs + (q) * (1 - filterprobs)) * (f(k) - norm_const) - level }
      }
      VaR = uniroot(g, interval = c(-1, 0))$root
      return(VaR)
    }
  }
}

# Function to calculate the Value at Risk (VaR) for a three-regime MSCMGARCH model
VaR_MSCMGARCH_3 <- function(portfolio_weights, H, level, type1, type2, par, filterprobs) {
  p1t = filterprobs[1]
  p2t = filterprobs[2]
  p3t = filterprobs[3]
  correlation_matrix_1 = cor_mat(par[7:9], type1)
  correlation_matrix_2 = cor_mat(par[10:12], type2)
  p11 = par[1]
  p12 = par[2]
  p13 = 1 - p11 - p12
  p21 = par[3]
  p22 = par[4]
  p23 = 1 - p21 - p22
  p32 = par[5]
  p33 = par[6]
  p31 = 1 - p33 - p32
  SD = sqrt(portfolio_weights %*% H %*% portfolio_weights)
  portfolio_1 = portfolio_weights %*% eigen_value_decomposition(H) %*% solve(t(chol(correlation_matrix_1)))
  portfolio_2 = portfolio_weights %*% eigen_value_decomposition(H) %*% solve(t(chol(correlation_matrix_2)))
  copula_par1 = par[7:9]
  copula_par2 = par[10:12]
  f = function(z) {
    VaR = pnorm(z, sd = SD)
    return(VaR)
  }
  f1 = function(k) {
    type = type1
    if (sum(portfolio_1 == 0) > 1) {
      if (sum(portfolio_1) > 0) {
        VaR = pnorm(k, sd = sum(portfolio_1))
      } else {
        VaR = 1 - pnorm(k, sd = sum(portfolio_1))
      }
      return(VaR)
    }
    if (portfolio_1[3] == 0) {
      if (portfolio_1[1] > 0) {
        return(cubintegrate(function(u) { copulaH2(pnorm(k / portfolio_1[1] - u * portfolio_1[2] / portfolio_1[1]), pnorm(u), copula_par1[1], type1[1]) * dnorm(u) }, lower = -2, upper = 2)$integral)
      } else {
        return((1 - cubintegrate(function(u) { copulaH2(pnorm(k / portfolio_1[1] - u * portfolio_1[2] / portfolio_1[1]), pnorm(u), copula_par1[1], type1[1]) * dnorm(u) }, lower = -2, upper = 2)$integral))
      }
    }
    if (portfolio_1[2] == 0) {
      if (portfolio[1] > 0) {
        return(cubintegrate(function(u) { return(copulaH2(pnorm(k / portfolio_1[1] - u * portfolio_1[3] / portfolio_1[1]), pnorm(u), copula_par1[2], type1[2]) * dnorm(u)) }, lower = -2, upper = 2)$integral)
      } else {
        return((1 - cubintegrate(function(u) { return(copulaH2(pnorm(k / portfolio_1[1] - u * portfolio_1[3] / portfolio_1[1]), pnorm(u), copula_par1[2], type1[2]) * dnorm(u)) }, lower = -2, upper = 2)$integral))
      }
    }
    if (portfolio_1[1] == 0) {
      if (portfolio_1[2] > 0) {
        return(cubintegrate(function(u) { copulaH2(pnorm(k / portfolio_1[2] - u * portfolio_1[3] / portfolio_1[2]), pnorm(u), copula_par1[3], type1[3]) * dnorm(u) }, lower = -2, upper = 2)$integral)
      } else {
        return((1 - cubintegrate(function(u) { copulaH2(pnorm(k / portfolio_1[2] - u * portfolio_1[3] / portfolio_1[2]), pnorm(u), copula_par1[3], type1[3]) * dnorm(u) }, lower = -2, upper = 2)$integral))
      }
    } else {
      if (portfolio_1[1] > 0) {
        f = function(k) { integral(function(x, y) { VineH2(pnorm(k / portfolio_1[1] - x * portfolio_1[2] / portfolio_1[1] - y * portfolio_1[3] / portfolio_1[1]), pnorm(x), pnorm(y), copula_par1, type1) * dnorm(x) * dnorm(y) }, bounds = list(x = c(-3, 3), y = c(-3, 3)), absTol = 1e-3)$value }
        norm_const = f(-1)
        return(f(k) - norm_const)
      } else {
        f = function(k) { 1 - integral(function(x, y) { VineH2(pnorm(k / portfolio_1[1] - x * portfolio_1[2] / portfolio_1[1] - y * portfolio_1[3] / portfolio_1[1]), pnorm(x), pnorm(y), copula_par1, type1) * dnorm(x) * dnorm(y) }, bounds = list(x = c(-3, 3), y = c(-3, 3)), absTol = 1e-3)$value }
        norm_const = f(-1)
        return(f(k) - norm_const)
      }
    }
  }
  f2 = function(k) {
    type = type2
    if (sum(portfolio_2 == 0) > 1) {
      if (sum(portfolio_2) > 0) {
        VaR = pnorm(k, sd = sum(portfolio_2))
      } else {
        VaR = 1 - pnorm(k, sd = sum(portfolio_2))
      }
      return(VaR)
    }
    if (portfolio_2[3] == 0) {
      if (portfolio_2[1] > 0) {
        return(cubintegrate(function(u) { copulaH2(pnorm(k / portfolio_2[1] - u * portfolio_2[2] / portfolio_2[1]), pnorm(u), copula_par2[1], type2[1]) * dnorm(u) }, lower = -2, upper = 2)$integral)
      } else {
        return((1 - cubintegrate(function(u) { copulaH2(pnorm(k / portfolio_2[1] - u * portfolio_2[2] / portfolio_2[1]), pnorm(u), copula_par2[1], type2[1]) * dnorm(u) }, lower = -2, upper = 2)$integral))
      }
    }
    if (portfolio_2[2] == 0) {
      if (portfolio[1] > 0) {
        return(cubintegrate(function(u) { return(copulaH2(pnorm(k / portfolio_2[1] - u * portfolio_2[3] / portfolio_2[1]), pnorm(u), copula_par2[2], type2[2]) * dnorm(u)) }, lower = -2, upper = 2)$integral)
      } else {
        return((1 - cubintegrate(function(u) { return(copulaH2(pnorm(k / portfolio_2[1] - u * portfolio_2[3] / portfolio_2[1]), pnorm(u), copula_par2[2], type2[2]) * dnorm(u)) }, lower = -2, upper = 2)$integral))
      }
    }
    if (portfolio_2[1] == 0) {
      if (portfolio_2[2] > 0) {
        return(cubintegrate(function(u) { copulaH2(pnorm(k / portfolio_2[2] - u * portfolio_2[3] / portfolio_2[2]), pnorm(u), copula_par2[3], type2[3]) * dnorm(u) }, lower = -2, upper = 2)$integral)
      } else {
        return((1 - cubintegrate(function(u) { copulaH2(pnorm(k / portfolio_2[2] - u * portfolio_2[3] / portfolio_2[2]), pnorm(u), copula_par2[3], type2[3]) * dnorm(u) }, lower = -2, upper = 2)$integral))
      }
    } else {
      if (portfolio_2[1] > 0) {
        f = function(k) { integral(function(x, y) { VineH2(pnorm(k / portfolio_2[1] - x * portfolio_2[2] / portfolio_2[1] - y * portfolio_2[3] / portfolio_2[1]), pnorm(x), pnorm(y), copula_par2, type2) * dnorm(x) * dnorm(y) }, bounds = list(x = c(-3, 3), y = c(-3, 3)), absTol = 1e-3)$value }
        norm_const = f(-1)
        return(f(k) - norm_const)
      } else {
        f = function(k) { 1 - integral(function(x, y) { VineH2(pnorm(k / portfolio_2[1] - x * portfolio_2[2] / portfolio_2[1] - y * portfolio_2[3] / portfolio_2[1]), pnorm(x), pnorm(y), copula_par2, type2) * dnorm(x) * dnorm(y) }, bounds = list(x = c(-3, 3), y = c(-3, 3)), absTol = 1e-3)$value }
        norm_const = f(-1)
        return(f(k) - norm_const)
      }
    }
  }
  final = function(k) {
    return((p1t * p11 + p2t * p21 + p3t * p31) * f(k) + (p1t * p12 + p2t * p22 + p3t * p32) * f1(k) + (p1t * p13 + p2t * p23 + p3t * p33) * f2(k) - level)
  }
  VaR = uniroot(final, interval = c(-1, 0))$root
  return(VaR)
}

# Function to calculate the Expected Shortfall (ES) for a three-regime MSCMGARCH model
ES_MSCMGARCH_3 <- function(portfolio_weights, H, level, type1, type2, par, filterprobs, VaR) {
  p1t = filterprobs[1]
  p2t = filterprobs[2]
  p3t = filterprobs[3]
  correlation_matrix_1 = cor_mat(par[7:9], type1)
  correlation_matrix_2 = cor_mat(par[10:12], type2)
  p11 = par[1]
  p12 = par[2]
  p13 = 1 - p11 - p12
  p21 = par[3]
  p22 = par[4]
  p23 = 1 - p21 - p22
  p32 = par[5]
  p33 = par[6]
  p31 = 1 - p33 - p32
  SD = sqrt(portfolio_weights %*% H %*% portfolio_weights)
  portfolio_1 = portfolio_weights %*% eigen_value_decomposition(H) %*% solve(t(chol(correlation_matrix_1)))
  portfolio_2 = portfolio_weights %*% eigen_value_decomposition(H) %*% solve(t(chol(correlation_matrix_2)))
  copula_par1 = par[7:9]
  copula_par2 = par[10:12]
  f3 = function(k) {
    VaR = integrate(function(z) { z * dnorm(z, sd = SD) }, lower = -Inf, upper = k, rel.tol = 1e-15)$value
    return(VaR)
  }
  f1 = function(k) {
    if (sum(portfolio_1 == 0) > 1) {
      f = function(x) {
        if (sum(portfolio_1) > 0) {
          VaR = pnorm(x, sd = sum(portfolio_1))
        } else {
          VaR = 1 - pnorm(x, sd = sum(portfolio_1))
        }
        return(VaR)
      }
      return(-integral(f, bounds = list(x = c(-Inf, k)), absTol = 1e-3)$value / f(k))
    }
    if (portfolio_1[3] == 0) {
      f = function(x) {
        if (portfolio_1[1] > 0) {
          return(cubintegrate(function(u) { copulaH2(pnorm(x / portfolio_1[1] - u * portfolio_1[2] / portfolio_1[1]), pnorm(u), copula_par1[1], type1[1]) * dnorm(u) }, lower = -2, upper = 2)$integral)
        } else {
          return((1 - cubintegrate(function(u) { copulaH2(pnorm(x / portfolio_1[1] - u * portfolio_1[2] / portfolio_1[1]), pnorm(u), copula_par1[1], type1[1]) * dnorm(u) }, lower = -2, upper = 2)$integral))
        }
      }
      return(-integral(f, bounds = list(x = c(-1, k)), absTol = 1e-3)$value / f(k))
    }
    if (portfolio_1[2] == 0) {
      f = function(x) {
        if (portfolio[1] > 0) {
          return(cubintegrate(function(u) { return(copulaH2(pnorm(x / portfolio_1[1] - u * portfolio_1[3] / portfolio_1[1]), pnorm(u), copula_par1[2], type1[2]) * dnorm(u)) }, lower = -2, upper = 2)$integral)
        } else {
          return((1 - cubintegrate(function(u) { return(copulaH2(pnorm(x / portfolio_1[1] - u * portfolio_1[3] / portfolio_1[1]), pnorm(u), copula_par1[2], type1[2]) * dnorm(u)) }, lower = -2, upper = 2)$integral))
        }
      }
      return(-integral(f, bounds = list(x = c(-1, k)), absTol = 1e-3)$value / f(k))
    }
    if (portfolio_1[1] == 0) {
      f = function(x) {
        if (portfolio_1[2] > 0) {
          return(cubintegrate(function(u) { copulaH2(pnorm(x / portfolio_1[2] - u * portfolio_1[3] / portfolio_1[2]), pnorm(u), copula_par1[3], type1[3]) * dnorm(u) }, lower = -2, upper = 2)$integral)
        } else {
          return((1 - cubintegrate(function(u) { copulaH2(pnorm(x / portfolio_1[2] - u * portfolio_1[3] / portfolio_1[2]), pnorm(u), copula_par1[3], type1[3]) * dnorm(u) }, lower = -2, upper = 2)$integral))
        }
      }
      return(-integral(f, bounds = list(x = c(-1, k)), absTol = 1e-3)$value / f(k))
    } else {
      f = function(x) { integral(function(y, z) { VineCopula(c(pnorm(x / portfolio_1[1] - y * portfolio_1[2] / portfolio_1[1] - z * portfolio_1[3] / portfolio_1[1]), pnorm(y), pnorm(z)), copula_par1, type1) * dnorm(x / portfolio_1[1] - y * portfolio_1[2] / portfolio_1[1] - z * portfolio_1[3] / portfolio_1[1]) * dnorm(y) * dnorm(z) }, bounds = list(y = c(-1, 1), z = c(-1, 1)), absTol = 1e-3)$value }
      norm_const = integral(function(z) { f(z) }, bounds = list(z = c(-1, k)), absTol = 1e-3)$value
      ES = integral(function(z) { z * (f(z)) }, bounds = list(z = c(-1, k)), absTol = 1e-3)$value / norm_const
      return(ES)
    }
  }
  f2 = function(k) {
    portfolio_1 = portfolio_2
    type1 = type2
    copula_par1 = copula_par2
    if (sum(portfolio_1 == 0) > 1) {
      f = function(x) {
        if (sum(portfolio_1) > 0) {
          VaR = pnorm(x, sd = sum(portfolio_1))
        } else {
          VaR = 1 - pnorm(x, sd = sum(portfolio_1))
        }
        return(VaR)
      }
      return(-integral(f, bounds = list(x = c(-Inf, k)), absTol = 1e-3)$value / f(k))
    }
    if (portfolio_1[3] == 0) {
      f = function(x) {
        if (portfolio_1[1] > 0) {
          return(cubintegrate(function(u) { copulaH2(pnorm(x / portfolio_1[1] - u * portfolio_1[2] / portfolio_1[1]), pnorm(u), copula_par1[1], type1[1]) * dnorm(u) }, lower = -2, upper = 2)$integral)
        } else {
          return((1 - cubintegrate(function(u) { copulaH2(pnorm(x / portfolio_1[1] - u * portfolio_1[2] / portfolio_1[1]), pnorm(u), copula_par1[1], type1[1]) * dnorm(u) }, lower = -2, upper = 2)$integral))
        }
      }
      return(-integral(f, bounds = list(x = c(-1, k)), absTol = 1e-3)$value / f(k))
    }
    if (portfolio_1[2] == 0) {
      f = function(x) {
        if (portfolio[1] > 0) {
          return(cubintegrate(function(u) { return(copulaH2(pnorm(x / portfolio_1[1] - u * portfolio_1[3] / portfolio_1[1]), pnorm(u), copula_par1[2], type1[2]) * dnorm(u)) }, lower = -2, upper = 2)$integral)
        } else {
          return((1 - cubintegrate(function(u) { return(copulaH2(pnorm(x / portfolio_1[1] - u * portfolio_1[3] / portfolio_1[1]), pnorm(u), copula_par1[2], type1[2]) * dnorm(u)) }, lower = -2, upper = 2)$integral))
        }
      }
      return(-integral(f, bounds = list(x = c(-1, k)), absTol = 1e-3)$value / f(k))
    }
    if (portfolio_1[1] == 0) {
      f = function(x) {
        if (portfolio_1[2] > 0) {
          return(cubintegrate(function(u) { copulaH2(pnorm(x / portfolio_1[2] - u * portfolio_1[3] / portfolio_1[2]), pnorm(u), copula_par1[3], type1[3]) * dnorm(u) }, lower = -2, upper = 2)$integral)
        } else {
          return((1 - cubintegrate(function(u) { copulaH2(pnorm(x / portfolio_1[2] - u * portfolio_1[3] / portfolio_1[2]), pnorm(u), copula_par1[3], type1[3]) * dnorm(u) }, lower = -2, upper = 2)$integral))
        }
      }
      return(-integral(f, bounds = list(x = c(-1, k)), absTol = 1e-3)$value / f(k))
    } else {
      f = function(x) { integral(function(y, z) { VineCopula(c(pnorm(x / portfolio_1[1] - y * portfolio_1[2] / portfolio_1[1] - z * portfolio_1[3] / portfolio_1[1]), pnorm(y), pnorm(z)), copula_par1, type1) * dnorm(x / portfolio_1[1] - y * portfolio_1[2] / portfolio_1[1] - z * portfolio_1[3] / portfolio_1[1]) * dnorm(y) * dnorm(z) }, bounds = list(y = c(-1, 1), z = c(-1, 1)), absTol = 1e-3)$value }
      norm_const = integral(function(z) { f(z) }, bounds = list(z = c(-1, k)), absTol = 1e-3)$value
      ES = integral(function(z) { z * (f(z)) }, bounds = list(z = c(-1, k)), absTol = 1e-3)$value / norm_const
      return(ES)
    }
  }
  return((p1t * p11 + p2t * p21 + p3t * p31) * f3(VaR) / pnorm(VaR, sd = SD) + (p1t * p12 + p2t * p22 + p3t * p32) * f1(VaR) + (p1t * p13 + p2t * p23 + p3t * p33) * f2(VaR))
}

# Function to calculate the u statistic for a three-regime MSCMGARCH model
u_MSCMGARCH_3 <- function(portfolio_weights, H, level, type1, type2, par, filterprobs, VaR, portfolio_return) {
  p1t = filterprobs[1]
  p2t = filterprobs[2]
  p3t = filterprobs[3]
  correlation_matrix_1 = cor_mat(par[7:9], type1)
  correlation_matrix_2 = cor_mat(par[10:12], type2)
  p11 = par[1]
  p12 = par[2]
  p13 = 1 - p11 - p12
  p21 = par[3]
  p22 = par[4]
  p23 = 1 - p21 - p22
  p32 = par[5]
  p33 = par[6]
  p31 = 1 - p33 - p32
  SD = sqrt(portfolio_weights %*% H %*% portfolio_weights)
  portfolio_1 = portfolio_weights %*% eigen_value_decomposition(H) %*% solve(t(chol(correlation_matrix_1)))
  portfolio_2 = portfolio_weights %*% eigen_value_decomposition(H) %*% solve(t(chol(correlation_matrix_2)))
  copula_par1 = par[7:9]
  copula_par2 = par[10:12]
  f = function(z) {
    VaR = pnorm(z, sd = SD)
    return(VaR)
  }
  f1 = function(k) {
    type = type1
    if (sum(portfolio_1 == 0) > 1) {
      if (sum(portfolio_1) > 0) {
        VaR = pnorm(k, sd = sum(portfolio_1))
      } else {
        VaR = 1 - pnorm(k, sd = sum(portfolio_1))
      }
      return(VaR)
    }
    if (portfolio_1[3] == 0) {
      if (portfolio_1[1] > 0) {
        return(cubintegrate(function(u) { copulaH2(pnorm(k / portfolio_1[1] - u * portfolio_1[2] / portfolio_1[1]), pnorm(u), copula_par1[1], type1[1]) * dnorm(u) }, lower = -2, upper = 2)$integral)
      } else {
        return((1 - cubintegrate(function(u) { copulaH2(pnorm(k / portfolio_1[1] - u * portfolio_1[2] / portfolio_1[1]), pnorm(u), copula_par1[1], type1[1]) * dnorm(u) }, lower = -2, upper = 2)$integral))
      }
    }
    if (portfolio_1[2] == 0) {
      if (portfolio[1] > 0) {
        return(cubintegrate(function(u) { return(copulaH2(pnorm(k / portfolio_1[1] - u * portfolio_1[3] / portfolio_1[1]), pnorm(u), copula_par1[2], type1[2]) * dnorm(u)) }, lower = -2, upper = 2)$integral)
      } else {
        return((1 - cubintegrate(function(u) { return(copulaH2(pnorm(k / portfolio_1[1] - u * portfolio_1[3] / portfolio_1[1]), pnorm(u), copula_par1[2], type1[2]) * dnorm(u)) }, lower = -2, upper = 2)$integral))
      }
    }
    if (portfolio_1[1] == 0) {
      if (portfolio_1[2] > 0) {
        return(cubintegrate(function(u) { copulaH2(pnorm(k / portfolio_1[2] - u * portfolio_1[3] / portfolio_1[2]), pnorm(u), copula_par1[3], type1[3]) * dnorm(u) }, lower = -2, upper = 2)$integral)
      } else {
        return((1 - cubintegrate(function(u) { copulaH2(pnorm(k / portfolio_1[2] - u * portfolio_1[3] / portfolio_1[2]), pnorm(u), copula_par1[3], type1[3]) * dnorm(u) }, lower = -2, upper = 2)$integral))
      }
    } else {
      if (portfolio_1[1] > 0) {
        f = function(k) { integral(function(x, y) { VineH2(pnorm(k / portfolio_1[1] - x * portfolio_1[2] / portfolio_1[1] - y * portfolio_1[3] / portfolio_1[1]), pnorm(x), pnorm(y), copula_par1, type1) * dnorm(x) * dnorm(y) }, bounds = list(x = c(-3, 3), y = c(-3, 3)))$value }
        norm_const = f(-1)
        g = return((integral(function(x, y) { VineH2(pnorm(k / portfolio_1[1] - x * portfolio_1[2] / portfolio_1[1] - y * portfolio_1[3] / portfolio_1[1]), pnorm(x), pnorm(y), copula_par1, type1) * dnorm(x) * dnorm(y) }, bounds = list(x = c(-3, 3), y = c(-3, 3)))$value - norm_const))
      } else {
        f = function(k) { 1 - integral(function(x, y) { VineH2(pnorm(k / portfolio_1[1] - x * portfolio_1[2] / portfolio_1[1] - y * portfolio_1[3] / portfolio_1[1]), pnorm(x), pnorm(y), copula_par1, type1) * dnorm(x) * dnorm(y) }, bounds = list(x = c(-3, 3), y = c(-3, 3)))$value }
        norm_const = f(-1)
        g = return((1 - integral(function(x, y) { VineH2(pnorm(k / portfolio_1[1] - x * portfolio_1[2] / portfolio_1[1] - y * portfolio_1[3] / portfolio_1[1]), pnorm(x), pnorm(y), copula_par1, type1) * dnorm(x) * dnorm(y) }, bounds = list(x = c(-3, 3), y = c(-3, 3)))$value - norm_const))
      }
      return(g)
    }
  }
  f2 = function(k) {
    type = type2
    if (sum(portfolio_2 == 0) > 1) {
      if (sum(portfolio_2) > 0) {
        VaR = pnorm(k, sd = sum(portfolio_2))
      } else {
        VaR = 1 - pnorm(k, sd = sum(portfolio_2))
      }
      return(VaR)
    }
    if (portfolio_2[3] == 0) {
      if (portfolio_2[1] > 0) {
        return(cubintegrate(function(u) { copulaH2(pnorm(k / portfolio_2[1] - u * portfolio_2[2] / portfolio_2[1]), pnorm(u), copula_par2[1], type2[1]) * dnorm(u) }, lower = -2, upper = 2)$integral)
      } else {
        return((1 - cubintegrate(function(u) { copulaH2(pnorm(k / portfolio_2[1] - u * portfolio_2[2] / portfolio_2[1]), pnorm(u), copula_par2[1], type2[1]) * dnorm(u) }, lower = -2, upper = 2)$integral))
      }
    }
    if (portfolio_2[2] == 0) {
      if (portfolio[1] > 0) {
        return(cubintegrate(function(u) { return(copulaH2(pnorm(k / portfolio_2[1] - u * portfolio_2[3] / portfolio_2[1]), pnorm(u), copula_par2[2], type2[2]) * dnorm(u)) }, lower = -2, upper = 2)$integral)
      } else {
        return((1 - cubintegrate(function(u) { return(copulaH2(pnorm(k / portfolio_2[1] - u * portfolio_2[3] / portfolio_2[1]), pnorm(u), copula_par2[2], type2[2]) * dnorm(u)) }, lower = -2, upper = 2)$integral))
      }
    }
    if (portfolio_2[1] == 0) {
      if (portfolio_2[2] > 0) {
        return(cubintegrate(function(u) { copulaH2(pnorm(k / portfolio_2[2] - u * portfolio_2[3] / portfolio_2[2]), pnorm(u), copula_par2[3], type2[3]) * dnorm(u) }, lower = -2, upper = 2)$integral)
      } else {
        return((1 - cubintegrate(function(u) { copulaH2(pnorm(k / portfolio_2[2] - u * portfolio_2[3] / portfolio_2[2]), pnorm(u), copula_par2[3], type2[3]) * dnorm(u) }, lower = -2, upper = 2)$integral))
      }
    } else {
      if (portfolio_2[1] > 0) {
        f = function(k) { integral(function(x, y) { VineH2(pnorm(k / portfolio_2[1] - x * portfolio_2[2] / portfolio_2[1] - y * portfolio_2[3] / portfolio_2[1]), pnorm(x), pnorm(y), copula_par2, type2) * dnorm(x) * dnorm(y) }, bounds = list(x = c(-3, 3), y = c(-3, 3)), absTol = 1e-3)$value }
        norm_const = f(-1)
        g = return((integral(function(x, y) { VineH2(pnorm(k / portfolio_2[1] - x * portfolio_2[2] / portfolio_2[1] - y * portfolio_2[3] / portfolio_2[1]), pnorm(x), pnorm(y), copula_par2, type2) * dnorm(x) * dnorm(y) }, bounds = list(x = c(-3, 3), y = c(-3, 3)), absTol = 1e-3)$value - norm_const))
      } else {
        f = function(k) { 1 - integral(function(x, y) { VineH2(pnorm(k / portfolio_2[1] - x * portfolio_2[2] / portfolio_2[1] - y * portfolio_2[3] / portfolio_2[1]), pnorm(x), pnorm(y), copula_par2, type2) * dnorm(x) * dnorm(y) }, bounds = list(x = c(-3, 3), y = c(-3, 3)), absTol = 1e-3)$value }
        norm_const = f(-1)
        return((1 - integral(function(x, y) { VineH2(pnorm(k / portfolio_2[1] - x * portfolio_2[2] / portfolio_2[1] - y * portfolio_2[3] / portfolio_2[1]), pnorm(x), pnorm(y), copula_par2, type2) * dnorm(x) * dnorm(y) }, bounds = list(x = c(-3, 3), y = c(-3, 3)), absTol = 1e-3)$value - norm_const))
      }
    }
  }
  final = function(k) {
    return((p1t * p11 + p2t * p21 + p3t * p31) * f(k) + (p1t * p12 + p2t * p22 + p3t * p32) * f1(k) + (p1t * p13 + p2t * p23 + p3t * p33) * f2(k))
  }
  if (portfolio_return <= VaR) {
    return(final(portfolio_return))
  }
  return(0)
}
