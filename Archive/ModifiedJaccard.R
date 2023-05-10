### The issue - Jaccard is symmetric such that the distance between com1 & com2 does not depend on the direction i.e. com1 -> com2 or com2 -> com1
# This matters because in the ocean there might be a difference in the difficulty to go from com1->com2 or vice versa due to currents etc.

# What about a modified jaccard that measures overlap form the perspective of com1->com2 ? 

#Some communities to play with
com1 <- c(1,1,1,1,0,0,0,0,1,1,1,1,0,0,0,0,0,0,0,0)
com2 <- c(0,0,0,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1)
com2 <- c(1,1,1,1,0,0,0,0)
com2 <- c(0,0,1,1,1,1,1,1)
com1 <- c(0,0,0,0,0,1,1,1)
com2 <- c(1,1,1,1,0,0,0,0)

#Jaccard
sum(com1&com2)/sum(com1|com2)
#Modified Jaccard
sum(com1&com2)/sum(com1)


#The modified distance between two sites considers... 
## what proportion of spp in com1 are also found in com2
## 1 if all spp in com1 are in com2 (different to vanilla jaccard)
## 0 if no spp in com1 are in com2 (same as vanilla jaccard)