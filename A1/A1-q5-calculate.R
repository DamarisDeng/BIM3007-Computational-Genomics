s <- 542
a0 <- 0.4309
b0 <- 0.5691
a <- 48*a0
a
b <- 48*b0
b
AB <- 134+a
Ab <- 154+b
aB <- 82+b
ab <- 76+a
print(c(AB, Ab, aB, ab))
print(c(AB/s, Ab/s, aB/s, ab/s))
P1 <- AB*ab/(AB*ab+Ab*aB)
P2 <- Ab*aB/(AB*ab+Ab*aB)
P1
P2
