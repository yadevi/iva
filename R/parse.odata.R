
parse.odata <- function(oformula, odata){

  if(!("Formula" %in% class(oformula))){
    if("formula" %in% class(oformula)){
      oformula<-Formula(oformula)
    }else{
      stop("oformula should be of class \"formula\"")
    }
  }

  if(!('id' %in% colnames(odata))){
    stop('Cannot find a column \'id\' in odata')
  }

  if(any(duplicated(odata$id))){
    stop('ID in odata is not unique')
  }

  rownames(odata) <- odata$id

  mf <- model.frame(oformula, na.action = na.pass, data = odata, rhs=1:2, lhs=1, drop=FALSE)
  ox <- model.matrix(oformula, mf, rhs=1, drop=F)[, -1, drop = FALSE]
  og <- model.matrix(oformula, mf, rhs=2, drop=F)[, -1, drop = FALSE]

  od <- model.part(oformula, mf, lhs=1, drop=FALSE)

  var.d <- colnames(od)
  var.ox <- colnames(ox)
  var.g <- colnames(og)
  odata <- data.frame(od, ox, og,
                      stringsAsFactors = FALSE,
                      check.names = FALSE)

  list(odata = odata, var.d = var.d, var.ox = var.ox, var.g = var.g)

}
