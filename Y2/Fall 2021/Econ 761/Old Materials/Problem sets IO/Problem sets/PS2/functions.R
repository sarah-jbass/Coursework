lerner=function(p, mc){
  return((p-mc)/p)
}

HHI=function(q,Q){
  si=q/Q
  return(sum(si^2))
}