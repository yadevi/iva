
parse.edata <- function(eformula, edata){

  if(!("Formula" %in% class(eformula))){
    if("formula" %in% class(eformula)){
      eformula<-Formula(eformula)
    }else{
      stop("rformula should be of class \"formula\"")
    }
  }

  if(!('id' %in% colnames(edata))){
    stop('Cannot find a column \'id\' in edata')
  }

  if(any(duplicated(edata$id))){
    stop('ID in edata is not unique')
  }

  rownames(edata) <- edata$id

  mf <- model.frame(eformula, na.action = na.pass, data = edata, rhs=1:2, lhs=1:2, drop=FALSE)
  ed <- model.part(eformula, mf, lhs=2, drop=F)[, , drop = FALSE]
  ex <- model.matrix(eformula, mf, rhs=1, drop=F)[, -1, drop = FALSE]
  eg <- model.matrix(eformula, mf, rhs=2, drop=F)[, -1, drop = FALSE]
  z <- model.part(eformula, mf, lhs=1, drop=FALSE)[, , drop = FALSE]

  var.d <- colnames(ed)
  var.z <- colnames(z)
  var.ex <- colnames(ex)
  var.g <- colnames(eg)

  edata <- data.frame(z, ed, ex, eg,
                      stringsAsFactors = FALSE,
                      check.names = FALSE)

  list(edata = edata, var.z = var.z, var.ex = var.ex, var.g = var.g, var.d = var.d)

}
