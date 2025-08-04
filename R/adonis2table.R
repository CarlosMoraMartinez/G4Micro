adonis2table <- function(mod1, perm_disp=NULL, ancap=NULL, var=""){
  aux <- data.frame(
    variable = ifelse(var == "", rownames(mod1)[1], var),
    DF_var=mod1$Df[1],
    DF_Residual = mod1$Df[2],
    DF_Total = mod1$Df[3],
    SumOfSQs_var = mod1$SumOfSqs[1],
    SumOfSQs_Residual =mod1$SumOfSqs[2],
    SumOfSQs_Total = mod1$SumOfSqs[3],
    R2_var = mod1$R2[1],
    R2_Residual = mod1$R2[2],
    R2_Total = mod1$R2[3],
    F_statistic = mod1$F[1],
    P = mod1$`Pr(>F)`[1]
  )
  if(!is.null(perm_disp)){
    aux$perm_disp_P = perm_disp$tab[1, "Pr(>F)"]
    aux$perm_disp_F = perm_disp$tab[1, "F"]
    aux$perm_disp_npermuts = perm_disp$tab[1, "N.Perm"]
  }else{
    aux$perm_disp_P = NA
    aux$perm_disp_F = NA
    aux$perm_disp_npermuts = NA
  }

  if(!is.null(ancap)){
    aux$capscaleanova_P = ancap[1, "Pr(>F)"]
    aux$capscaleanopva_F = ancap[1, "F"]
  }else{
    aux$capscaleanova_P = NA
    aux$capscaleanopva_F = NA
  }

  return(aux)
}
