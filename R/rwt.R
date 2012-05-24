rwt <-
function(d){
  
  tRank <- function(x){
    n<-length(x)
    rt <- rank(x, ties.method= "random")
  }
  
  rd <- apply(t(d),1,tRank)
  return(rd)
}
