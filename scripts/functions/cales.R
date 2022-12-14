# From https://github.com/wt12318/NeoEnrichment/
cales <- function(mutation_dt,neoantigen_list,type){
  a <- sum(mutation_dt[mutation_dt$index %in% neoantigen_list,"rank"])
  if(type=="I"){
    b <- (nrow(mutation_dt)-length(neoantigen_list))
    mutation_dt$re <- ifelse(mutation_dt$index %in% neoantigen_list,((mutation_dt$rank)/a),-(1/b))
    
  }else{
    b <- sum(mutation_dt[!(mutation_dt$index %in% neoantigen_list),"rank"])
    mutation_dt$re <- ifelse(mutation_dt$index %in% neoantigen_list,((mutation_dt$rank)/a),-((mutation_dt$rank)/b))
  }
  mutation_dt$cum_re <- cumsum(mutation_dt$re)
  es <- max(0,mutation_dt$cum_re)-abs(min(mutation_dt$cum_re,0))
}
