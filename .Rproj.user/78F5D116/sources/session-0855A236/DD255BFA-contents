
#' Title
#'
#' @param phi0 probability
#' @param phi probability
#' @param la scale parameter
#' @param th discrete parameter
#'
#' @return requires a value
#' @export
#'
#' @examples
#' CZIGPFx(0.5,0.5,1,0.1)
CZIGPFx = function(phi0, phi, la, th)
{
  m = length(phi)
  df = matrix(0,3*m+1,3*m+1)
  a0 = prod(phi+(1-phi)*exp(-la))
  a = rep(0,m)
  d = matrix(0,m,m)
  for(i in 1:m)
  {
    a[i] = prod( phi[-i]+(1-phi[-i])*exp(-la[-i]))
    for(j in 1:m)
    {
      d[j,i] = d[i,j] = prod( phi[-c(i,j)]+(1-phi[-c(i,j)])*exp(-la[-c(i,j)]) )
    }
  }
  r1 = phi0+(1-phi0)*a0
  df[1,1] = -(1-a0)^2/r1 - (1-a0)/(1-phi0)
  df[2:(m+1),1] = df[1,2:(m+1)] = -(1-exp(-la))*a/r1
  df[(m+2):(2*m+1),1] = df[1,(m+2):(2*m+1)] = (1-phi)*exp(-la)*a + (1-phi0)*(1-a0)*(1-phi)*exp(-la)*a/r1
  df[(2*m+2):(3*m+1),1] = df[1,(2*m+2):(3*m+1)] = 0
  for(i in 1:m)
  {
    for(j in 1:m)
    {
      df[j+1,i+1] = df[i+1,j+1] = (1-phi0)*(1-exp(-la[i]))*(1-exp(-la[j]))*d[i,j] - (1-phi0)^2*(1-exp(-la[i]))*(1-exp(-la[j]))*a[i]*a[j]/r1
      df[j+m+1,i+1] = df[i+1,j+m+1] = (1-phi0)^2*(1-exp(-la[i]))*a[i]*(1-phi[j])*exp(-la[j])*a[j]/r1 - (1-phi0)*(1-exp(-la[i]))*(1-phi[j])*exp(-la[j])*d[i,j]
      df[j+2*m+1,i+1] = df[i+1,j+2*m+1] =0
      df[j+m+1,i+m+1] = df[i+m+1,j+m+1] = (1-phi0)*(1-phi[i])*(1-phi[j])*exp(-la[i])*exp(-la[j])*d[i,j]-(1-phi0)^2*(1-phi[i])*(1-phi[j])*exp(-la[i])*exp(-la[j])*a[i]*a[j]/r1
      df[j+2*m+1,i+m+1] = df[i+m+1,j+2*m+1] = 0
      df[j+2*m+1,i+2*m+1] = df[i+2*m+1,j+2*m+1] =0
    }

    df[i+1,i+1] = -(1-phi0)^2*(1-exp(-la[i]))^2*a[i]^2/r1-(1-r1)*(1-exp(-la[i]))^2/(phi[i]+(1-phi[i])*exp(-la[i]))-(1-r1)*(1-exp(-la[i]))/(1-phi[i])
    df[i+m+1,i+m+1] = (1-phi0)*(1-phi[i])*exp(-la[i])*a[i] - (1-phi0)^2*(1-phi[i])^2*exp(-2*la[i])*a[i]^2/r1 + phi[i]*(1-phi[i])*exp(-la[i])*(1-r1)/(phi[i]+(1-phi[i])*exp(-la[i]))+(1-phi[i])*(1-exp(-la[i]))*(1-r1)*(th[i]/(la[i]+2*th[i])-1/la[i])
    df[i+2*m+1,i+2*m+1] = -(1-phi[i])*(1-exp(-la[i]))*(1-r1)*( 2*la[i]/(la[i]+2*th[i])+la[i]/(1-th[i]) )

    df[i+m+1,i+1] = df[i+1,i+m+1] = (1-phi0)*exp(-la[i])*a[i] + (1-phi0)^2*(1-exp(-la[i]))*(1-phi[i])*exp(-la[i])*a[i]^2/r1 + (1-r1)*exp(-la[i])/(phi[i]+(1-phi[i])*exp(-la[i]))
    df[i+2*m+1,i+1] = df[i+1,i+2*m+1] = 0
    df[i+2*m+1,i+m+1] = df[i+m+1,i+2*m+1] = -(1-phi[i])*(1-exp(-la[i]))*(1-r1)*la[i]/(la[i]+2*th[i])
  }

  return(df)
}
