load("C:\\Users\\Finlay\\Downloads\\VL.RData")

attach(VL)
actual_diff = diff(by(titer,type,median))

compare_med = function(x,factor,permutations){
  x_temp = x
  storage = vector(length=permutations)
  for(i in 1:permutations){
    median_temp = tapply(sample(x_temp),factor,median)
    storage[i] = median_temp[1]-median_temp[2]
  }
  return(as.randtest(sim=storage,actual_diff))
}

compare_med(titer,type,1000)