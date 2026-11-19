# CARE1 1.1.1: public functions as.record() and estN().
# Original estimation formulas retained; helper functions are internal.

.care1_estimates <- function(z, se = FALSE, nboot = 200) {
  z <- as.vector(unlist(z, use.names = FALSE))
  if (!is.numeric(z) || length(z) != 7L && length(z) != 15L && length(z) != 31L ||
      anyNA(z) || any(!is.finite(z)) || any(z < 0) || any(z != floor(z))) {
    stop("z must be a nonnegative integer count vector of length 7, 15, or 31.",
         call. = FALSE)
  }
  if (!is.logical(se) || length(se) != 1L || is.na(se)) {
    stop("se must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.numeric(nboot) || length(nboot) != 1L || is.na(nboot) ||
      !is.finite(nboot) || nboot < 2 || nboot != floor(nboot)) {
    stop("nboot must be an integer greater than one.", call. = FALSE)
  }
  
  k <- match(length(z), c(7L, 15L, 31L)) + 2L
  pat <- vapply(seq_len(2^k - 1L), function(v) {
    as.integer(intToBits(v))[k:1L]
  }, integer(k))
  n_list <- as.numeric(pat %*% z)
  if (any(n_list == 0)) {
    stop("Each list must contain at least one case.", call. = FALSE)
  }
  M <- sum(z)
  
  funs <- list("Nhat-0" = .care1_indep,
               "Nhat-1" = .care1_lsc,
               "Nhat" = .care1_hsc)
  point <- vapply(funs, function(f) as.numeric(f(z)), numeric(1L))
  result <- matrix(NA_real_, nrow = 3L, ncol = 4L,
                   dimnames = list(names(funs), c("est", "se", "cil", "ciu")))
  result[, "est"] <- point
  
  if (se) {
    for (i in seq_along(funs)) {
      nhat <- point[i]
      if (!is.finite(nhat) || nhat < M || nhat <= 0 ||
          nhat > .Machine$integer.max) {
        warning("Bootstrap cannot be performed for ", names(funs)[i],
                ": estimated population size is invalid or below M.",
                call. = FALSE)
        next
      }
      zstar <- .care1_bootstrap(z, nhat, nboot)
      b <- apply(zstar, 2L, function(x) {
        tryCatch(
          suppressWarnings(as.numeric(funs[[i]](x))),
          error = function(e) NA_real_
        )
      })
      if (any(!is.finite(b))) {
        warning("Non-finite bootstrap estimates for ", names(funs)[i],
                "; uncertainty is reported as NA.", call. = FALSE)
        next
      }
      s <- stats::sd(b)
      result[i, "se"] <- s
      if (nhat > M && is.finite(s)) {
        factor <- exp(1.96 * sqrt(log1p(s^2 / (nhat - M)^2)))
        result[i, "cil"] <- M + (nhat - M) / factor
        result[i, "ciu"] <- M + (nhat - M) * factor
      } else if (nhat == M && s == 0) {
        result[i, c("cil", "ciu")] <- M
      } else {
        warning("Log-transformed CI is undefined for ", names(funs)[i],
                " when the estimate equals M and bootstrap SE > 0.",
                call. = FALSE)
      }
    }
  }
  
  return(result)
}

# Equivalent to sampling Nhat individuals from observed categories and
# the estimated unseen category. Only observed categories are returned.
.care1_bootstrap <- function(z, nhat, nboot) {
  size <- as.integer(round(nhat))
  missing <- size - sum(z)
  if (missing < 0) stop("Rounded population estimate is smaller than M.",
                        call. = FALSE)
  prob <- c(missing, z) / size
  stats::rmultinom(nboot, size = size, prob = prob)[-1L, , drop = FALSE]
}


# Existing CARE1 formula for Nhat-0 (kept unchanged).
.care1_indep <- function(z){
  if(length(z)==7){
    z001=z[1];	z010=z[2];	z011=z[3];	z100=z[4];	z101=z[5];	z110=z[6];	z111=z[7]
    n1=z100+z101+z110+z111
    n2=z010+z011+z110+z111
    n3=z001+z011+z101+z111
    chat=1-(1/3)*(z100/n1+z010/n2+z001/n3)
    d=(1/3)*(z100+z101+z010+z011+z110+z111
             +z100+z110+z001+z011+z101+z111
             +z010+z110+z001+z101+z011+z111)
    nhat0=d/chat
    nhat0
  }
  else if(length(z)==15){
    z0001=z[1];	z0010=z[2]; z0011=z[3]; z0100=z[4]; z0101=z[5];
    z0110=z[6]; z0111=z[7]; z1000=z[8]; z1001=z[9]; z1010=z[10];
    z1011=z[11]; z1100=z[12]; z1101=z[13]; z1110=z[14]; z1111=z[15]
    
    n1=z1000+z1001+z1010+z1100+z1011+z1110+z1101+z1111
    n2=z0100+z0101+z0110+z0111+z1100+z1110+z1101+z1111
    n3=z0010+z0011+z0110+z0111+z1010+z1011+z1110+z1111
    n4=z0001+z0011+z0101+z0111+z1001+z1011+z1101+z1111
    chat=1-(1/4)*(z1000/n1+z0100/n2+z0010/n3+z0001/n4)
    d=(1/4)*(z1000+z1001+z1010+z1011+z1100+z1101+z0100+z0101+z0110+z0111+z0010+z0011+z1110+z1111  
             +z1000+z1010+z1001+z1011+z1100+z1110+z0001+z0011+z0101+z0111+z0100+z0110+z1101+z1111
             +z1000+z1100+z1001+z1101+z1010+z1110+z0001+z0101+z0010+z0110+z0011+z0111+z1011+z1111
             +z0001+z1001+z0010+z1010+z0011+z1011+z0100+z1100+z0101+z1101+z0110+z1110+z0111+z1111)
    nhat0=d/chat
    nhat0
  }
  else if(length(z)==31){
    z00001=z[1]; z00010=z[2]; z00011=z[3]; z00100=z[4];
    z00101=z[5]; z00110=z[6]; z00111=z[7]; z01000=z[8]; 
    z01001=z[9]; z01010=z[10]; z01011=z[11]; z01100=z[12];
    z01101=z[13]; z01110=z[14]; z01111=z[15]; z10000=z[16]; 
    z10001=z[17]; z10010=z[18]; z10011=z[19]; z10100=z[20];
    z10101=z[21]; z10110=z[22]; z10111=z[23]; z11000=z[24];
    z11001=z[25]; z11010=z[26]; z11011=z[27]; z11100=z[28];
    z11101=z[29]; z11110=z[30]; z11111=z[31]
    
    n1=z10000+z10001+z10010+z10011+z10100+z10101+z10110+z10111+z11000+z11001+z11010+z11011+z11100+z11101+z11110+z11111
    n2=z01000+z01001+z01010+z01011+z01100+z01101+z01110+z01111+z11000+z11001+z11010+z11011+z11100+z11101+z11110+z11111
    n3=z00100+z00101+z00110+z00111+z01100+z01101+z01110+z01111+z10100+z10101+z10110+z10111+z11100+z11101+z11110+z11111
    n4=z00010+z00011+z00110+z00111+z01010+z01011+z01110+z01111+z10010+z10011+z10110+z10111+z11010+z11011+z11110+z11111
    n5=z00001+z00011+z00101+z00111+z01001+z01011+z01101+z01111+z10001+z10011+z10101+z10111+z11001+z11011+z11101+z11111
    chat=1-(1/5)*(z10000/n1+z01000/n2+z00100/n3+z00010/n4+z00001/n5)
    
    d=(1/5)*(z00010+z00011+z00100+z00101+z00110+z00111+z01000+z01001+z01010+z01011+z01100+z01101
             +z01110+z01111+z10000+z10001+z10010+z10011+z10100+z10101+z10110+z10111+z11000+z11001
             +z11010+z11011+z11100+z11101+z11110+z11111
             +z00001+z00011+z00100+z00110+z00101+z00111+z01000+z01010+z01001+z01011+z01100+z01110
             +z01101+z01111+z10000+z10010+z10001+z10011+z10100+z10110+z10101+z10111+z11000+z11010
             +z11001+z11011+z11100+z11110+z11101+z11111
             +z00001+z00101+z00010+z00110+z00011+z00111+z01000+z01100+z01001+z01101+z01010+z01110
             +z01011+z01111+z10000+z10100+z10001+z10101+z10010+z10110+z10011+z10111+z11000+z11100
             +z11001+z11101+z11010+z11110+z11011+z11111
             +z00001+z01001+z00010+z01010+z00011+z01011+z00100+z01100+z00101+z01101+z00110+z01110
             +z00111+z01111+z10000+z11000+z10001+z11001+z10010+z11010+z10011+z11011+z10100+z11100
             +z10101+z11101+z10110+z11110+z10111+z11111
             +z00001+z10001+z00010+z10010+z00011+z10011+z00100+z10100+z00101+z10101+z00110+z10110
             +z00111+z10111+z01000+z11000+z01001+z11001+z01010+z11010+z01011+z11011+z01100+z11100
             +z01101+z11101+z01110+z11110+z01111+z11111)
    nhat0=d/chat
    nhat0
  }
  else
    "Invalid argument error! Please input valid ascertainment records."
}



# Existing CARE1 formula for Nhat-1 (kept unchanged).
.care1_lsc <- function(z){
  if(length(z)==7){
    z001=z[1];	z010=z[2];	z011=z[3];	z100=z[4];	z101=z[5];	z110=z[6];	z111=z[7]
    n1=z100+z101+z110+z111
    n2=z010+z011+z110+z111
    n3=z001+z011+z101+z111
    chat=1-(1/3)*(z100/n1+z010/n2+z001/n3)
    d=(1/3)*(z100+z101+z010+z011+z110+z111
             +z100+z110+z001+z011+z101+z111
             +z010+z110+z001+z101+z011+z111)
    n0=d/chat
    n=n0+((z100+z110+z010+z110)*(n0*(z110+z111)/(n1*n2)-1)
          +(z100+z101+z001+z101)*(n0*(z101+z111)/(n1*n3)-1)
          +(z010+z011+z001+z011)*(n0*(z011+z111)/(n2*n3)-1))/(3*chat)
    r12=n*(z110+z111)/(n1*n2)-1
    r13=n*(z101+z111)/(n1*n3)-1
    r23=n*(z011+z111)/(n2*n3)-1
    y=((z100+z110+z010+z110)*r12+(z100+z101+z001+z101)*r13+(z010+z011+z001+z011)*r23)
    nhat1=d/chat+1/(3*chat)*y
    nhat1
  }
  else if(length(z)==15){
    z0001=z[1];	z0010=z[2]; z0011=z[3]; z0100=z[4]; z0101=z[5];
    z0110=z[6]; z0111=z[7]; z1000=z[8]; z1001=z[9]; z1010=z[10];
    z1011=z[11]; z1100=z[12]; z1101=z[13]; z1110=z[14]; z1111=z[15]
    
    n1=z1000+z1001+z1010+z1100+z1011+z1110+z1101+z1111
    n2=z0100+z0101+z0110+z0111+z1100+z1110+z1101+z1111
    n3=z0010+z0011+z0110+z0111+z1010+z1011+z1110+z1111
    n4=z0001+z0011+z0101+z0111+z1001+z1011+z1101+z1111
    chat=1-(1/4)*(z1000/n1+z0100/n2+z0010/n3+z0001/n4)
    d=(1/4)*(z1000+z1001+z1010+z1011+z1100+z1101+z0100+z0101+z0110+z0111+z0010+z0011+z1110+z1111  
             +z1000+z1010+z1001+z1011+z1100+z1110+z0001+z0011+z0101+z0111+z0100+z0110+z1101+z1111
             +z1000+z1100+z1001+z1101+z1010+z1110+z0001+z0101+z0010+z0110+z0011+z0111+z1011+z1111
             +z0001+z1001+z0010+z1010+z0011+z1011+z0100+z1100+z0101+z1101+z0110+z1110+z0111+z1111)			
    
    n0=d/chat
    n=n0+((z1000+z1100+z0100+z1100)*(n0*(z1100+z1101+z1110+z1111)/(n1*n2)-1)
          +(z1000+z1010+z0010+z1010)*(n0*(z1010+z1011+z1110+z1111)/(n1*n3)-1)
          +(z1000+z1001+z0001+z1001)*(n0*(z1001+z1011+z1101+z1111)/(n1*n4)-1)
          +(z0100+z0110+z0010+z0110)*(n0*(z0110+z0111+z1110+z1111)/(n2*n3)-1)
          +(z0100+z0101+z0001+z0101)*(n0*(z0101+z0111+z1101+z1111)/(n2*n4)-1)
          +(z0010+z0011+z0001+z0011)*(n0*(z0011+z0111+z1011+z1111)/(n3*n4)-1))/(4*chat)
    
    r12=n*(z1100+z1101+z1110+z1111)/(n1*n2)-1                                          
    r13=n*(z1010+z1011+z1110+z1111)/(n1*n3)-1
    r14=n*(z1001+z1011+z1101+z1111)/(n1*n4)-1
    r23=n*(z0110+z0111+z1110+z1111)/(n2*n3)-1
    r24=n*(z0101+z0111+z1101+z1111)/(n2*n4)-1
    r34=n*(z0011+z0111+z1011+z1111)/(n3*n4)-1
    y=((z1000+z1100+z0100+z1100)*r12+(z1000+z1010+z0010+z1010)*r13+(z1000+z1001+z0001+z1001)*r14
       +(z0100+z0110+z0010+z0110)*r23+(z0100+z0101+z0001+z0101)*r24+(z0010+z0011+z0001+z0011)*r34)
    nhat1=d/chat+1/(4*chat)*y
    nhat1
  }
  else if(length(z)==31){
    z00001=z[1]; z00010=z[2]; z00011=z[3]; z00100=z[4];
    z00101=z[5]; z00110=z[6]; z00111=z[7]; z01000=z[8]; 
    z01001=z[9]; z01010=z[10]; z01011=z[11]; z01100=z[12];
    z01101=z[13]; z01110=z[14]; z01111=z[15]; z10000=z[16]; 
    z10001=z[17]; z10010=z[18]; z10011=z[19]; z10100=z[20];
    z10101=z[21]; z10110=z[22]; z10111=z[23]; z11000=z[24];
    z11001=z[25]; z11010=z[26]; z11011=z[27]; z11100=z[28];
    z11101=z[29]; z11110=z[30]; z11111=z[31]
    
    n1=z10000+z10001+z10010+z10011+z10100+z10101+z10110+z10111+z11000+z11001+z11010+z11011+z11100+z11101+z11110+z11111
    n2=z01000+z01001+z01010+z01011+z01100+z01101+z01110+z01111+z11000+z11001+z11010+z11011+z11100+z11101+z11110+z11111
    n3=z00100+z00101+z00110+z00111+z01100+z01101+z01110+z01111+z10100+z10101+z10110+z10111+z11100+z11101+z11110+z11111
    n4=z00010+z00011+z00110+z00111+z01010+z01011+z01110+z01111+z10010+z10011+z10110+z10111+z11010+z11011+z11110+z11111
    n5=z00001+z00011+z00101+z00111+z01001+z01011+z01101+z01111+z10001+z10011+z10101+z10111+z11001+z11011+z11101+z11111
    chat=1-(1/5)*(z10000/n1+z01000/n2+z00100/n3+z00010/n4+z00001/n5)
    
    d=(1/5)*(z00010+z00011+z00100+z00101+z00110+z00111+z01000+z01001+z01010+z01011+z01100+z01101
             +z01110+z01111+z10000+z10001+z10010+z10011+z10100+z10101+z10110+z10111+z11000+z11001
             +z11010+z11011+z11100+z11101+z11110+z11111
             +z00001+z00011+z00100+z00110+z00101+z00111+z01000+z01010+z01001+z01011+z01100+z01110
             +z01101+z01111+z10000+z10010+z10001+z10011+z10100+z10110+z10101+z10111+z11000+z11010
             +z11001+z11011+z11100+z11110+z11101+z11111
             +z00001+z00101+z00010+z00110+z00011+z00111+z01000+z01100+z01001+z01101+z01010+z01110
             +z01011+z01111+z10000+z10100+z10001+z10101+z10010+z10110+z10011+z10111+z11000+z11100
             +z11001+z11101+z11010+z11110+z11011+z11111
             +z00001+z01001+z00010+z01010+z00011+z01011+z00100+z01100+z00101+z01101+z00110+z01110
             +z00111+z01111+z10000+z11000+z10001+z11001+z10010+z11010+z10011+z11011+z10100+z11100
             +z10101+z11101+z10110+z11110+z10111+z11111
             +z00001+z10001+z00010+z10010+z00011+z10011+z00100+z10100+z00101+z10101+z00110+z10110
             +z00111+z10111+z01000+z11000+z01001+z11001+z01010+z11010+z01011+z11011+z01100+z11100
             +z01101+z11101+z01110+z11110+z01111+z11111)
    
    n0=d/chat
    n=n0+((z10000+z11000+z01000+z11000)*(n0*(z11000+z11001+z11010+z11011+z11100+z11101+z11110+z11111)/(n1*n2)-1)
          +(z10000+z10100+z00100+z10100)*(n0*(z10100+z10101+z10110+z10111+z11100+z11101+z11110+z11111)/(n1*n3)-1)
          +(z10000+z10010+z00010+z10010)*(n0*(z10010+z10011+z10110+z10111+z11011+z11010+z11110+z11111)/(n1*n4)-1)
          +(z10000+z10001+z00001+z10001)*(n0*(z10001+z10011+z10101+z10111+z11001+z11011+z11101+z11111)/(n1*n5)-1)
          +(z01000+z01100+z00100+z01100)*(n0*(z01100+z01101+z01110+z01111+z11100+z11101+z11110+z11111)/(n2*n3)-1)
          +(z01000+z01010+z00010+z01010)*(n0*(z01010+z01011+z01110+z01111+z11010+z11011+z11110+z11111)/(n2*n4)-1)
          +(z01000+z01001+z00001+z01001)*(n0*(z01001+z01011+z01101+z01111+z11001+z11011+z11101+z11111)/(n2*n5)-1)
          +(z00100+z00110+z00010+z00110)*(n0*(z00110+z00111+z01110+z01111+z10110+z10111+z11110+z11111)/(n3*n4)-1)
          +(z00100+z00101+z00001+z00101)*(n0*(z00101+z00111+z01101+z01111+z10101+z10111+z11101+z11111)/(n3*n5)-1)
          +(z00010+z00011+z00001+z00011)*(n0*(z00011+z00111+z01011+z01111+z10011+z10111+z11011+z11111)/(n4*n5)-1))/(5*chat)
    
    r12=n*(z11000+z11001+z11010+z11011+z11100+z11101+z11110+z11111)/(n1*n2)-1                                        
    r13=n*(z10100+z10101+z10110+z10111+z11100+z11101+z11110+z11111)/(n1*n3)-1
    r14=n*(z10010+z10011+z10110+z10111+z11011+z11010+z11110+z11111)/(n1*n4)-1
    r15=n*(z10001+z10011+z10101+z10111+z11001+z11011+z11101+z11111)/(n1*n5)-1
    r23=n*(z01100+z01101+z01110+z01111+z11100+z11101+z11110+z11111)/(n2*n3)-1
    r24=n*(z01010+z01011+z01110+z01111+z11010+z11011+z11110+z11111)/(n2*n4)-1
    r25=n*(z01001+z01011+z01101+z01111+z11001+z11011+z11101+z11111)/(n2*n5)-1
    r34=n*(z00110+z00111+z01110+z01111+z10110+z10111+z11110+z11111)/(n3*n4)-1
    r35=n*(z00101+z00111+z01101+z01111+z10101+z10111+z11101+z11111)/(n3*n5)-1
    r45=n*(z00011+z00111+z01011+z01111+z10011+z10111+z11011+z11111)/(n4*n5)-1
    
    y=((z10000+z11000+z01000+z11000)*r12+(z10000+z10100+z00100+z10100)*r13
       +(z10000+z10010+z00010+z10010)*r14+(z10000+z10001+z00001+z10001)*r15
       +(z01000+z01100+z00100+z01100)*r23+(z01000+z01010+z00010+z01010)*r24
       +(z01000+z01001+z00001+z01001)*r25+(z00100+z00110+z00010+z00110)*r34
       +(z00100+z00101+z00001+z00101)*r35+(z00010+z00011+z00001+z00011)*r45)
    
    nhat1=d/chat+1/(5*chat)*y
    nhat1
  }
  else
    "Invalid argument error! Please input valid ascertainment records."
}

# Existing CARE1 formula for Nhat (kept unchanged).
.care1_hsc <- function(z){
  if(length(z)==7){
    z001=z[1];	z010=z[2];	z011=z[3];	z100=z[4];	z101=z[5];	z110=z[6];	z111=z[7]
    n1=z100+z101+z110+z111
    n2=z010+z011+z110+z111
    n3=z001+z011+z101+z111
    chat=1-(1/3)*(z100/n1+z010/n2+z001/n3)
    M=z001+z010+z011+z100+z101+z110+z111
    x=3*chat-((z100+z110+z010+z110)*(z110+z111)/(n1*n2)
              +(z100+z101+z001+z101)*(z101+z111)/(n1*n3)
              +(z001+z011+z010+z011)*(z011+z111)/(n2*n3))
    y=(z011+z111+z101+z111+z110+z111)/x
    nhat=ifelse(x<0.0001 || y<M, M, y)
    nhat
  }
  else if(length(z)==15){
    z0001=z[1];	z0010=z[2]; z0011=z[3]; z0100=z[4]; z0101=z[5];
    z0110=z[6]; z0111=z[7]; z1000=z[8]; z1001=z[9]; z1010=z[10];
    z1011=z[11]; z1100=z[12]; z1101=z[13]; z1110=z[14]; z1111=z[15]
    
    n1=z1000+z1001+z1010+z1100+z1011+z1110+z1101+z1111
    n2=z0100+z0101+z0110+z0111+z1100+z1110+z1101+z1111
    n3=z0010+z0011+z0110+z0111+z1010+z1011+z1110+z1111
    n4=z0001+z0011+z0101+z0111+z1001+z1011+z1101+z1111
    chat=1-(1/4)*(z1000/n1+z0100/n2+z0010/n3+z0001/n4)
    M=z0001+z0010+z0011+z0100+z0101+z0110+z0111+z1000+z1001+z1010+z1100+z1011+z1110+z1101+z1111
    x=4*chat-((z1000+z1100+z0100+z1100)*(z1100+z1101+z1110+z1111)/(n1*n2) 
              +(z1000+z1010+z0010+z1010)*(z1010+z1011+z1110+z1111)/(n1*n3)
              +(z1000+z1001+z0001+z1001)*(z1001+z1011+z1101+z1111)/(n1*n4)
              +(z0100+z0110+z0010+z0110)*(z0110+z0111+z1110+z1111)/(n2*n3)
              +(z0100+z0101+z0001+z0101)*(z0101+z0111+z1101+z1111)/(n2*n4)
              +(z0010+z0011+z0001+z0011)*(z0011+z0111+z1011+z1111)/(n3*n4))
    y=(z1110+z1111+z1101+z1111+z1011+z1111+z0111+z1111+z1100+z1101+z1010+z1011+z0110+z0111+z1100+z1110
       +z1001+z1011+z0101+z0111+z1010+z1110+z1001+z1101+z0011+z0111+z0110+z1110+z0101+z1101+z0011+z1011)/x
    nhat=ifelse(x<0.0001 || y<M, M, y)
    nhat
  }
  else if(length(z)==31){
    z00001=z[1]; z00010=z[2]; z00011=z[3]; z00100=z[4];
    z00101=z[5]; z00110=z[6]; z00111=z[7]; z01000=z[8]; 
    z01001=z[9]; z01010=z[10]; z01011=z[11]; z01100=z[12];
    z01101=z[13]; z01110=z[14]; z01111=z[15]; z10000=z[16]; 
    z10001=z[17]; z10010=z[18]; z10011=z[19]; z10100=z[20];
    z10101=z[21]; z10110=z[22]; z10111=z[23]; z11000=z[24];
    z11001=z[25]; z11010=z[26]; z11011=z[27]; z11100=z[28];
    z11101=z[29]; z11110=z[30]; z11111=z[31]
    
    n1=z10000+z10001+z10010+z10011+z10100+z10101+z10110+z10111+z11000+z11001+z11010+z11011+z11100+z11101+z11110+z11111
    n2=z01000+z01001+z01010+z01011+z01100+z01101+z01110+z01111+z11000+z11001+z11010+z11011+z11100+z11101+z11110+z11111
    n3=z00100+z00101+z00110+z00111+z01100+z01101+z01110+z01111+z10100+z10101+z10110+z10111+z11100+z11101+z11110+z11111
    n4=z00010+z00011+z00110+z00111+z01010+z01011+z01110+z01111+z10010+z10011+z10110+z10111+z11010+z11011+z11110+z11111
    n5=z00001+z00011+z00101+z00111+z01001+z01011+z01101+z01111+z10001+z10011+z10101+z10111+z11001+z11011+z11101+z11111
    chat=1-(1/5)*(z10000/n1+z01000/n2+z00100/n3+z00010/n4+z00001/n5)
    
    M=sum(z)
    
    x=5*chat-((z10000+z11000+z01000+z11000)*(z11000+z11001+z11010+z11011+z11100+z11101+z11110+z11111)/(n1*n2) 
              +(z10000+z10100+z00100+z10100)*(z10100+z10101+z10110+z10111+z11100+z11101+z11110+z11111)/(n1*n3)
              +(z10000+z10010+z00010+z10010)*(z10010+z10011+z10110+z10111+z11011+z11010+z11110+z11111)/(n1*n4)
              +(z10000+z10001+z00001+z10001)*(z10001+z10011+z10101+z10111+z11001+z11011+z11101+z11111)/(n1*n5)
              +(z01000+z01100+z00100+z01100)*(z01100+z01101+z01110+z01111+z11100+z11101+z11110+z11111)/(n2*n3)
              +(z01000+z01010+z00010+z01010)*(z01010+z01011+z01110+z01111+z11010+z11011+z11110+z11111)/(n2*n4)
              +(z01000+z01001+z00001+z01001)*(z01001+z01011+z01101+z01111+z11001+z11011+z11101+z11111)/(n2*n5)
              +(z00100+z00110+z00010+z00110)*(z00110+z00111+z01110+z01111+z10110+z10111+z11110+z11111)/(n3*n4)
              +(z00100+z00101+z00001+z00101)*(z00101+z00111+z01101+z01111+z10101+z10111+z11101+z11111)/(n3*n5)
              +(z00010+z00011+z00001+z00011)*(z00011+z00111+z01011+z01111+z10011+z10111+z11011+z11111)/(n4*n5))
    
    y=(z00110+z00111+z01010+z01011+z01100+z01101+z01110+z01111+z10010+z10011+z10100+z10101 
       +z10110+z10111+z11000+z11001+z11010+z11011+z11100+z11101+z11110+z11111
       +z00101+z00111+z01001+z01011+z01100+z01110+z01101+z01111+z10001+z10011+z10100+z10110
       +z10101+z10111+z11000+z11010+z11001+z11011+z11100+z11110+z11101+z11111
       +z00011+z00111+z01001+z01101+z01010+z01110+z01011+z01111+z10001+z10101+z10010+z10110
       +z10011+z10111+z11000+z11100+z11001+z11101+z11010+z11110+z11011+z11111
       +z00011+z01011+z00101+z01101+z00110+z01110+z00111+z01111+z10001+z11001+z10010+z11010
       +z10011+z11011+z10100+z11100+z10101+z11101+z10110+z11110+z10111+z11111
       +z00011+z10011+z00101+z10101+z00110+z10110+z00111+z10111+z01001+z11001+z01010+z11010
       +z01011+z11011+z01100+z11100+z01101+z11101+z01110+z11110+z01111+z11111)/x
    
    nhat=ifelse(x<0.0001 || y<M, M, y)
    nhat
  }
  else
    "Invalid argument error! Please input valid ascertainment records."
}



# Internal helper for 3, 4, or 5 ascertainment lists.
.care1_pattern_info <- function(z) {
  z <- as.vector(unlist(z, use.names = FALSE))
  if (!is.numeric(z) || !length(z) %in% c(7L, 15L, 31L) ||
      anyNA(z) || any(!is.finite(z)) || any(z < 0) ||
      any(z != floor(z))) {
    stop("z must contain 7, 15, or 31 nonnegative integer counts.", call. = FALSE)
  }
  k <- match(length(z), c(7L, 15L, 31L)) + 2L
  patterns <- t(vapply(seq_len(2^k - 1L), function(i) {
    as.integer(intToBits(i))[k:1L]
  }, integer(k)))
  colnames(patterns) <- paste0("n", seq_len(k))
  n <- as.numeric(crossprod(z, patterns))
  names(n) <- colnames(patterns)
  list(z = z, k = k, patterns = patterns, n = n, M = sum(z))
}

.care1_stat <- function(z) {
  a <- .care1_pattern_info(z)
  if (any(a$n == 0)) {
    stop("Every list must contain observations.", call. = FALSE)
  }
  singleton <- rowSums(a$patterns) == 1L
  D <- a$M - sum(a$z[singleton]) / a$k
  singleton_counts <- vapply(seq_len(a$k), function(i) {
    sum(a$z[singleton & a$patterns[, i] == 1L])
  }, numeric(1L))
  chat <- 1 - mean(singleton_counts / a$n)
  values <- c(a$n, M = a$M, D = round(D, 3), SC = round(chat, 3))
  matrix(as.numeric(values), nrow = 1L,
         dimnames = list(NULL, names(values)))
}

.care1_pair <- function(z) {
  a <- .care1_pattern_info(z)
  pairs <- utils::combn(seq_len(a$k), 2L)
  ans <- t(vapply(seq_len(ncol(pairs)), function(j) {
    i <- pairs[1L, j]; h <- pairs[2L, j]
    in_i <- a$patterns[, i] == 1L
    in_h <- a$patterns[, h] == 1L
    overlap <- sum(a$z[in_i & in_h])
    only_i <- sum(a$z[in_i & !in_h])
    only_h <- sum(a$z[!in_i & in_h])
    ni <- a$n[i]; nh <- a$n[h]
    Mij <- only_i + only_h + overlap
    pet <- if (overlap == 0) Inf else ni * nh / overlap
    chp <- (ni + 1) * (nh + 1) / (overlap + 1) - 1
    s <- sqrt((ni + 1) * (nh + 1) * only_i * only_h /
                ((overlap + 1)^2 * (overlap + 2)))
    if (chp > Mij && is.finite(s)) {
      factor <- exp(1.96 * sqrt(log1p(s^2 / (chp - Mij)^2)))
      cil <- Mij + (chp - Mij) / factor
      ciu <- Mij + (chp - Mij) * factor
    } else if (chp == Mij && s == 0) {
      cil <- ciu <- Mij
    } else {
      cil <- ciu <- NA_real_
    }
    c(Petersen = pet, Chapman = chp, se = s, cil = cil, ciu = ciu)
  }, numeric(5L)))
  rownames(ans) <- paste0("pa", pairs[1L, ], pairs[2L, ])
  ans
}

.care1_para <- function(z, nhat) {
  a <- .care1_pattern_info(z)
  if (!is.numeric(nhat) || length(nhat) != 1L || is.na(nhat) ||
      !is.finite(nhat) || nhat < 0) {
    stop("nhat must be a finite, nonnegative single population estimate.", call. = FALSE)
  }
  pairs <- utils::combn(seq_len(a$k), 2L)
  nm <- c(paste0("u", seq_len(a$k)),
          paste0("r", pairs[1L, ], pairs[2L, ]),
          if (a$k == 3L) "r123")
  if (nhat <= a$M) {
    out <- matrix(0, nrow = 1L, ncol = length(nm))
    colnames(out) <- nm
    return(out)
  }
  if (any(a$n == 0)) stop("Every list must contain observations.", call. = FALSE)
  u <- a$n / nhat
  r <- vapply(seq_len(ncol(pairs)), function(j) {
    i <- pairs[1L, j]; h <- pairs[2L, j]
    overlap <- sum(a$z[a$patterns[, i] == 1L & a$patterns[, h] == 1L])
    nhat * overlap / (a$n[i] * a$n[h]) - 1
  }, numeric(1L))
  if (a$k == 3L) {
    all_three <- sum(a$z[rowSums(a$patterns) == 3L])
    r123 <- nhat^2 * all_three / prod(a$n) - sum(r) - 1
    val <- c(u, r, r123)
  } else val <- c(u, r)
  out <- matrix(as.numeric(val), nrow = 1L)
  colnames(out) <- nm
  out
}


# Convert individual-level binary records to CARE1 categorical frequencies.
#' Transform Individual Data to CARE1 Input Format
#'
#' @description
#' Converts individual capture histories or presence/absence ascertainment
#' records into the categorical input format required by CARE1.
#'
#' @param x A matrix containing the observed individual capture histories
#'   or presence/absence ascertainment records.
#' @return A vector containing the number of individuals in the
#'   7, 15, or 31 nonzero capture-history categories.
#' @author T. C. Hsieh
#' @details
#' Each row of the input matrix represents an individual's capture history
#' or presence/absence ascertainment record across multiple sampling
#' occasions or lists. Each matrix element must be 0 (absence or
#' non-capture) or 1 (presence or capture).
#'
#' For example, consider the following matrix containing the
#' presence/absence records of 10 individuals across three lists. The
#' function classifies these 10 individuals into seven capture-history
#' categories: (001), (010), (011), (100), (101), (110), and (111),
#' and computes the number of individuals in each category. In this example,
#' three individuals were recorded in lists 1 and 3 but not in list 2
#' (category 101), while no individuals were recorded in all three lists
#' (category 111). Similar counts can be obtained for other categories.
#'
#' Data for 10 individuals:
#' \preformatted{
#'       List 1 List 2 List 3
#'  [1,]      0      1      0
#'  [2,]      1      0      1
#'  [3,]      0      1      0
#'  [4,]      0      1      1
#'  [5,]      1      0      1
#'  [6,]      1      0      1
#'  [7,]      0      0      1
#'  [8,]      0      1      0
#'  [9,]      0      1      0
#' [10,]      1      1      0
#' }
#'
#' The converted data are:
#' \preformatted{
#' 001 010 011 100 101 110 111
#'   1   4   1   0   3   1   0
#' }
#'
#' For four-list data, individual records are classified into 15
#' capture-history categories; for five-list data, they are classified
#' into 31 categories.
#' @examples
#' # Classify the capture histories of 100 individuals across three lists
#' # into seven categories.
#' set.seed(123)
#' x <- matrix(sample(0:1, 300, replace = TRUE), ncol = 3)
#' as.record(x)
#' @export
as.record <- function(x) {
  if (!is.matrix(x) && !is.data.frame(x)) {
    stop("x must be a matrix or data frame.", call. = FALSE)
  }
  x <- as.matrix(x)
  if (!is.numeric(x) && !is.logical(x)) {
    stop("x must contain only numeric/logical 0 and 1 values.", call. = FALSE)
  }
  if (!nrow(x) || !ncol(x) %in% 3:5 || anyNA(x) ||
      any(!(x %in% c(0, 1)))) {
    stop("x must have 3, 4, or 5 columns containing only 0 and 1.", call. = FALSE)
  }
  k <- ncol(x)
  values <- as.integer(x %*% 2^((k - 1L):0L))
  frequencies <- tabulate(values + 1L, nbins = 2^k)[-1L]
  pattern_names <- vapply(seq_len(2^k - 1L), function(v) {
    paste(rev(as.integer(intToBits(v))[seq_len(k)]), collapse = "")
  }, character(1L))
  stats::setNames(frequencies, pattern_names)
}

# Unified printed analysis; depends on the public functions defined above.

#' Population Size Estimation Using Capture-Recapture Models
#'
#' @description
#' Prints the number of identified cases in each list, optional pairwise
#' population size estimates, and three population size estimates using
#' the CARE1 sample coverage framework, together with parameter estimates.
#' @param z A vector in the categorical input format required by CARE1,
#'   or an output object returned by \code{as.record}.
#' @param se A logical variable indicating whether to compute bootstrap
#'   standard errors and log-transformed 95\% confidence intervals for
#'   the three population size estimates. Default is \code{TRUE}.
#'   This argument does not affect the analytical Chapman standard errors
#'   and confidence intervals.
#' @param nboot Number of bootstrap replications, an integer greater
#'   than one. Default is 100.
#' @param pair A logical variable indicating whether to compute and print
#'   estimates based on each pair of lists. Default is \code{TRUE}.
#' @return Invisibly returns a list containing:
#' \describe{
#'   \item{list_counts}{Number of identified cases in each list.}
#'   \item{pair_estimates}{Petersen and Chapman estimates, Chapman standard
#'     errors and log-transformed 95\% confidence intervals; NULL if pair is FALSE.}
#'   \item{population_estimates}{A matrix with rows Nhat-0, Nhat, Nhat-1
#'     and columns M, D, SC, est, se, cil, ciu.}
#'   \item{parameters}{Estimated mean ascertainment probabilities and
#'     coefficients of covariation for each population size estimate.}
#' }
#' M is the number of distinct observed individuals; D is the average
#' number of individuals in the union of all lists except one;
#' SC is estimated sample coverage.
#' \describe{
#'   \item{est}{Population size estimate.}
#'   \item{se}{Estimated standard error based on the bootstrap method.}
#'   \item{cil}{Lower limit of the 95\% confidence interval, calculated
#'     using a log-transformation.}
#'   \item{ciu}{Upper limit of the 95\% confidence interval, calculated
#'     using a log-transformation.}
#' }
#' When se is FALSE, the sample coverage uncertainty columns contain NA.
#' Definitions and an estimator selection guideline are also printed.
#'
#' @details
#' \describe{
#'   \item{Nhat-0}{Population size estimate under the assumption of
#'     independence among sources; see Eq. (9) of Chao et al. (2001).}
#'   
#'   \item{Nhat-1}{One-step population size estimate for cases with low
#'     sample coverage; see Eq. (17) of Chao et al. (2001), or Eq. (14)
#'     of Chao (2015) for an explicit formula for three-list data.}
#'   
#'   \item{Nhat}{Population size estimate for cases with sufficiently high
#'     sample coverage; see Eq. (16) of Chao et al. (2001), or Eq. (13)
#'     of Chao (2015) for an explicit formula for three-list data.}
#' }
#'
#'
#' @note
#' As a rough guideline, sample coverage is considered sufficiently high
#' when the estimated sample coverage (SC in the output) is at least 55\%
#' and the estimated bootstrap standard error of Nhat is less than one-third
#' of the population size estimate. Otherwise, Nhat-1 is recommended.
#'
#'
#' @references
#' Chao, A. (2015). Capture-recapture for human populations.
#' Wiley StatsRef: Statistics Reference Online. W. Piegorsch, Ed.
#' (John Wiley & Sons), 1-16.
#' \doi{10.1002/9781118445112.stat04855.pub2}
#'
#' Chao, A., Tsay, P. K., Lin, S. H., Shau, W. Y., and Chao, D. Y.
#' (2001). The applications of capture-recapture models to epidemiological
#' data. Statistics in Medicine, 20, 3123-3157.
#' \doi{10.1002/sim.996}
#'
#' @examples
#' data(Data_HAV)
#' estN(Data_HAV, se = FALSE)
#' estN(Data_HAV, se = FALSE, pair = FALSE)
#' \donttest{
#' set.seed(123)
#' estN(Data_HAV, se = TRUE, nboot = 100, pair = TRUE)
#' }
#' @export
estN <- function(z, se = TRUE, nboot = 100, pair = TRUE) {
  if (!is.logical(pair) || length(pair) != 1L || is.na(pair)) {
    stop("pair must be TRUE or FALSE.", call. = FALSE)
  }
  
  z <- as.vector(unlist(z, use.names = FALSE))
  a <- .care1_pattern_info(z)
  
  # 1. Compute population size estimates and summary statistics
  fit <- .care1_estimates(z, se = se, nboot = nboot)
  summary_stats <- .care1_stat(z)
  
  list_counts <- summary_stats[
    , paste0("n", seq_len(a$k)), drop = FALSE
  ]
  
  pairs <- utils::combn(seq_len(a$k), 2L)
  pair_estimates <- NULL
  if (pair) {
    pair_estimates <- .care1_pair(z)
    colnames(pair_estimates) <- c("Petersen", "Chapman", "se", "cil", "ciu")
    rownames(pair_estimates) <- paste0("pair(", pairs[1L, ], ",", pairs[2L, ], ")")
  }
  
  # 2. Preserve the original estimator display order
  kinds <- c("Nhat-0", "Nhat", "Nhat-1")
  
  pop_estimates <- cbind(
    M = rep(summary_stats[1L, "M"], length(kinds)),
    D = rep(summary_stats[1L, "D"], length(kinds)),
    SC = rep(summary_stats[1L, "SC"], length(kinds)),
    fit[kinds, , drop = FALSE]
  )
  
  rownames(pop_estimates) <- kinds
  colnames(pop_estimates) <- c(
    "M", "D", "SC", "est", "se", "cil", "ciu"
  )
  
  # 3. Parameter estimates
  n_parameters <- a$k + choose(a$k, 2L) +
    as.integer(a$k == 3L)
  
  parameters <- t(vapply(kinds, function(kind) {
    
    n_hat <- fit[kind, "est"]
    
    if (!is.finite(n_hat)) {
      return(rep(NA_real_, n_parameters))
    }
    
    as.numeric(.care1_para(z, n_hat))
    
  }, numeric(n_parameters)))
  
  colnames(parameters) <- c(
    paste0("u", seq_len(a$k)),
    paste0("r", pairs[1L, ], pairs[2L, ]),
    if (a$k == 3L) "r123"
  )
  
  rownames(parameters) <- kinds
  
  # ============================================================
  # PRINT RESULTS
  # ============================================================
  
  cat("Number of identified cases in each list:\n")
  print(as.data.frame(list_counts), row.names = FALSE)
  
  if (pair) {
    cat("\n ESTIMATES BASED ON EACH PAIR OF LISTS:\n")
    print(round(pair_estimates, 2))
    cat("\n    Note 1: Petersen and Chapman estimators provide population size estimates\n",
        "            based on pairs of lists. The Chapman standard errors are reported.\n",
        "    Note 2: A log-transformation is used to construct the 95% confidence\n",
        "            intervals for the Chapman estimates.\n", sep = "")
  }
  cat(if (pair) "\n SAMPLE COVERAGE APPROACH:\n" else "\n SAMPLE COVERAGE APPROACH:\n")
  print(round(pop_estimates, 2))
  if (!se) cat("Bootstrap standard errors and confidence intervals were not computed (se = FALSE).\n")
  
  # Diagnostic for the HSC estimator
  diag <- tryCatch(
    lee1(z),
    error = function(e) NA_real_
  )
  
  if (length(diag) == 1L &&
      is.finite(diag) &&
      diag <= 0) {
    
    cat(
      "\n",
      "    Warning: The HSC denominator is non-positive.\n",
      "             Nhat may be unstable. Consider using\n",
      "             Nhat-0 or Nhat-1.\n",
      sep = ""
    )
  }
  
  cat("\n")
  cat("    Parameter estimates:\n")
  print(round(parameters, 2))
  
  # ============================================================
  # COMPLETE DEFINITIONS
  # ============================================================
  
  cat(
    "\n",
    "DEFINITIONS for the sample coverage approach:\n",
    "\n",
    "M: Number of distinct individuals ascertained",
    " in at least one list.\n",
    "\n",
    "D: Average number of individuals in the union",
    " of all lists except one.\n",
    "\n",
    "SC: Estimated sample coverage.\n",
    "\n",
    "est: Population size estimate.\n",
    "\n",
    "se: Estimated standard error of the population",
    " size estimate based\n",
    "    on bootstrap replications.\n",
    "\n",
    "cil: Lower limit of the 95% confidence interval",
    " using a log-transformation.\n",
    "\n",
    "ciu: Upper limit of the 95% confidence interval",
    " using a log-transformation.\n",
    "\n",
    "Nhat-0: Population size estimate under the",
    " assumption of independence\n",
    "        among sources (Chao et al., 2001, Eq. 9).\n",
    "\n",
    "Nhat-1: One-step population size estimate for",
    " cases with low sample\n",
    "        coverage (Chao et al., 2001, Eq. 17).\n",
    "\n",
    "Nhat: Population size estimate for cases with",
    " sufficiently high\n",
    "      sample coverage (Chao et al., 2001, Eq. 16).\n",
    "\n",
    "u1, u2, ...: Estimated mean ascertainment",
    " probabilities depending\n",
    "             on the population size estimate.\n",
    "\n",
    "r12, r13, ...: Estimated coefficients of",
    " covariation (CCV)\n",
    "              depending on the population size estimate.\n",
    sep = ""
  )
  
  cat(
    "\n",
    "GUIDELINE FOR ESTIMATOR SELECTION:\n",
    "\n",
    "As a rough guideline, sample coverage is considered",
    " sufficiently\n",
    "high when the estimated sample coverage (SC in the output) is at least 55%\n",
    "and the estimated bootstrap",
    " standard error\n",
    "of Nhat is less than one-third of Nhat. Otherwise,",
    " Nhat-1 is\n",
    "recommended.\n",
    sep = ""
  )
  
  invisible(list(
    list_counts = list_counts,
    pair_estimates = pair_estimates,
    population_estimates = pop_estimates,
    parameters = parameters
  ))
}


# Original CARE1 denominator diagnostic; intentionally not exported.
lee1=function(z) {
  z=as.vector(unlist(z))
  lee3 <-
    function(z){
      z001=z[1];	z010=z[2];	z011=z[3];	z100=z[4];	z101=z[5];	z110=z[6];	z111=z[7]
      n1=z100+z101+z110+z111
      n2=z010+z011+z110+z111
      n3=z001+z011+z101+z111
      chat=1-(1/3)*(z100/n1+z010/n2+z001/n3)
      M=z001+z010+z011+z100+z101+z110+z111
      x=3*chat-((z100+z110+z010+z110)*(z110+z111)/(n1*n2)
                +(z100+z101+z001+z101)*(z101+z111)/(n1*n3)
                +(z001+z011+z010+z011)*(z011+z111)/(n2*n3))
      ifelse(x>0,1,0)
    }
  lee4 <-
    function(z){
      z0001=z[1];	z0010=z[2]; z0011=z[3]; z0100=z[4]; z0101=z[5];
      z0110=z[6]; z0111=z[7]; z1000=z[8]; z1001=z[9]; z1010=z[10];
      z1011=z[11]; z1100=z[12]; z1101=z[13]; z1110=z[14]; z1111=z[15]
      
      n1=z1000+z1001+z1010+z1100+z1011+z1110+z1101+z1111
      n2=z0100+z0101+z0110+z0111+z1100+z1110+z1101+z1111
      n3=z0010+z0011+z0110+z0111+z1010+z1011+z1110+z1111
      n4=z0001+z0011+z0101+z0111+z1001+z1011+z1101+z1111
      chat=1-(1/4)*(z1000/n1+z0100/n2+z0010/n3+z0001/n4)
      M=z0001+z0010+z0011+z0100+z0101+z0110+z0111+z1000+z1001+z1010+z1100+z1011+z1110+z1101+z1111
      x=4*chat-((z1000+z1100+z0100+z1100)*(z1100+z1101+z1110+z1111)/(n1*n2) 
                +(z1000+z1010+z0010+z1010)*(z1010+z1011+z1110+z1111)/(n1*n3)
                +(z1000+z1001+z0001+z1001)*(z1001+z1011+z1101+z1111)/(n1*n4)
                +(z0100+z0110+z0010+z0110)*(z0110+z0111+z1110+z1111)/(n2*n3)
                +(z0100+z0101+z0001+z0101)*(z0101+z0111+z1101+z1111)/(n2*n4)
                +(z0010+z0011+z0001+z0011)*(z0011+z0111+z1011+z1111)/(n3*n4))
      ifelse(x>0,1,0)
    }
  
  lee5 <-
    function(z){
      z00001=z[1]; z00010=z[2]; z00011=z[3]; z00100=z[4];
      z00101=z[5]; z00110=z[6]; z00111=z[7]; z01000=z[8]; 
      z01001=z[9]; z01010=z[10]; z01011=z[11]; z01100=z[12];
      z01101=z[13]; z01110=z[14]; z01111=z[15]; z10000=z[16]; 
      z10001=z[17]; z10010=z[18]; z10011=z[19]; z10100=z[20];
      z10101=z[21]; z10110=z[22]; z10111=z[23]; z11000=z[24];
      z11001=z[25]; z11010=z[26]; z11011=z[27]; z11100=z[28];
      z11101=z[29]; z11110=z[30]; z11111=z[31]
      
      n1=z10000+z10001+z10010+z10011+z10100+z10101+z10110+z10111+z11000+z11001+z11010+z11011+z11100+z11101+z11110+z11111
      n2=z01000+z01001+z01010+z01011+z01100+z01101+z01110+z01111+z11000+z11001+z11010+z11011+z11100+z11101+z11110+z11111
      n3=z00100+z00101+z00110+z00111+z01100+z01101+z01110+z01111+z10100+z10101+z10110+z10111+z11100+z11101+z11110+z11111
      n4=z00010+z00011+z00110+z00111+z01010+z01011+z01110+z01111+z10010+z10011+z10110+z10111+z11010+z11011+z11110+z11111
      n5=z00001+z00011+z00101+z00111+z01001+z01011+z01101+z01111+z10001+z10011+z10101+z10111+z11001+z11011+z11101+z11111
      chat=1-(1/5)*(z10000/n1+z01000/n2+z00100/n3+z00010/n4+z00001/n5)
      
      M=sum(z)
      
      x=5*chat-((z10000+z11000+z01000+z11000)*(z11000+z11001+z11010+z11011+z11100+z11101+z11110+z11111)/(n1*n2) 
                +(z10000+z10100+z00100+z10100)*(z10100+z10101+z10110+z10111+z11100+z11101+z11110+z11111)/(n1*n3)
                +(z10000+z10010+z00010+z10010)*(z10010+z10011+z10110+z10111+z11011+z11010+z11110+z11111)/(n1*n4)
                +(z10000+z10001+z00001+z10001)*(z10001+z10011+z10101+z10111+z11001+z11011+z11101+z11111)/(n1*n5)
                +(z01000+z01100+z00100+z01100)*(z01100+z01101+z01110+z01111+z11100+z11101+z11110+z11111)/(n2*n3)
                +(z01000+z01010+z00010+z01010)*(z01010+z01011+z01110+z01111+z11010+z11011+z11110+z11111)/(n2*n4)
                +(z01000+z01001+z00001+z01001)*(z01001+z01011+z01101+z01111+z11001+z11011+z11101+z11111)/(n2*n5)
                +(z00100+z00110+z00010+z00110)*(z00110+z00111+z01110+z01111+z10110+z10111+z11110+z11111)/(n3*n4)
                +(z00100+z00101+z00001+z00101)*(z00101+z00111+z01101+z01111+z10101+z10111+z11101+z11111)/(n3*n5)
                +(z00010+z00011+z00001+z00011)*(z00011+z00111+z01011+z01111+z10011+z10111+z11011+z11111)/(n4*n5))
      ifelse(x>0,1,0)
    }
  
  
  
  if(length(z)==7){lee3(z)}
  else if(length(z)==15){lee4(z)}
  else if(length(z)==31){lee5(z)}
}
