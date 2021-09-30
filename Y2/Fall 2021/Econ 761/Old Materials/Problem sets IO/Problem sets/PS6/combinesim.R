combinesim <- function(l1, l2, env = parent.frame()){
  nr <- l1$nr+l2$nr
  ph <- l1$ph+l2$ph
  return(list(nr = nr, ph = ph))
}