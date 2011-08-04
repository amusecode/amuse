      real function nfwdens(rad)

      include 'commonblocks'
 
      s = rad/a
      nfwdens = haloconst/(s**cusp)/((1.+s)**(outerhaloslope-cusp))
      return
      end
 
      real function nfwdensprime(rad)

      include 'commonblocks'
      real nfwdens

      s = rad/a
      nfwdensprime = -nfwdens(rad)/a*(outerhaloslope*s+cusp)/s/(1.+s)
      return
      end

      real function nfwdens2prime(rad)

      include 'commonblocks'
      real nfwdens

      s = rad/a
      nfwdens2prime = nfwdens(rad)/a/a*
     *     (cusp*(cusp+1.)+2.*(1.+outerhaloslope)*cusp*s+
     +     outerhaloslope*(1.+outerhaloslope)*s*s)/s/s/(1.+s)/(1.+s)
      return
      end

