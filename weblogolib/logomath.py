#!/usr/bin/env python

# -------------------------------- WebLogo --------------------------------

#  Copyright (c) 2003-2004 The Regents of the University of California.
#  Copyright (c) 2005 Gavin E. Crooks
#  Copyright (c) 2006-2011, The Regents of the University of California, through 
#  Lawrence Berkeley National Laboratory (subject to receipt of any required
#  approvals from the U.S. Dept. of Energy).  All rights reserved.

#  This software is distributed under the new BSD Open Source License.
#  <http://www.opensource.org/licenses/bsd-license.html>
#
#  Redistribution and use in source and binary forms, with or without 
#  modification, are permitted provided that the following conditions are met: 
#
#  (1) Redistributions of source code must retain the above copyright notice, 
#  this list of conditions and the following disclaimer. 
#
#  (2) Redistributions in binary form must reproduce the above copyright 
#  notice, this list of conditions and the following disclaimer in the 
#  documentation and or other materials provided with the distribution. 
#
#  (3) Neither the name of the University of California, Lawrence Berkeley 
#  National Laboratory, U.S. Dept. of Energy nor the names of its contributors 
#  may be used to endorse or promote products derived from this software 
#  without specific prior written permission. 
#
#  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
#  AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
#  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
#  ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE 
#  LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
#  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF 
#  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
#  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
#  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
#  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
#  POSSIBILITY OF SUCH DAMAGE. 

from math import log, sqrt, exp

from numpy import array, asarray, float64, ones, zeros, int32,all,any, shape
import numpy as na
from scipy import stats

from corebio.moremath import *

import random

# function: calculates p-values with binomial test
def binom_test(sample_cnt, composition):
    L = len(sample_cnt)
    S = sum(sample_cnt)
    all_cnt = [ S for i in range(L) ]
    probs = []
    for p in stats.binom.cdf(sample_cnt, all_cnt, composition):
        if p < 0.5: probs.append(-p)
        else:       probs.append(1 - p)
    return probs

# function: calculates p-values with t-test
def t_test(sample_cnt, composition):
    L = len(sample_cnt)
    S1 = asarray([sum(sample_cnt) for i in range(L)])
    S2 = asarray([sum(S1[0] * composition) for i in range(L)])
    c1 = sample_cnt
    c2 = S1 * composition
    m1 = c1 / S1
    m2 = c2 / S2
    v1 = c1 * (1 - m1) * (1 - m1) + (S1 - c1) * m1 * m1
    v2 = c2 * (1 - m2) * (1 - m2) + (S2 - c2) * m2 * m2

    df = S1 + S2 - asarray([2 for i in range(L)])
    var = (v1 + v2) / df
    probs = []
    for i in range(L):
         p = stats.t.cdf((m1[i] - m2[i]) / sqrt(var[i] * (1.0/S1[i] + 1.0/S2[i])), df[0])
         if p < 0.5: probs.append(-p)
         else:       probs.append(1 - p)
    return probs


# function: calculates p-values with z-score
def z_score(sample_cnt, composition):
    L = len(sample_cnt)
    S = sum(sample_cnt)
    bg_cnt = S * composition
    var = [ sqrt(i) for i in (ones(L, float64) - composition) * bg_cnt ]
    Z = (sample_cnt - bg_cnt)/var
    P = stats.norm.cdf(Z)
    sign = []
    for i in range(L):
        if Z[i] < 0:
            sign.append(-1);
        else       :
            P[i] = 1.0 - P[i]
            sign.append( 1)
    return list(P * sign)

class Dirichlet(object) :
    """The Dirichlet probability distribution. The Dirichlet is a continuous 
    multivariate probability distribution across non-negative unit length
    vectors. In other words, the Dirichlet is a probability distribution of 
    probability distributions. It is conjugate to the multinomial
    distribution and is widely used in Bayesian statistics.
    
    The Dirichlet probability distribution of order K-1 is 

     p(theta_1,...,theta_K) d theta_1 ... d theta_K = 
        (1/Z) prod_i=1,K theta_i^{alpha_i - 1} delta(1 -sum_i=1,K theta_i)

    The normalization factor Z can be expressed in terms of gamma functions:

      Z = {prod_i=1,K Gamma(alpha_i)} / {Gamma( sum_i=1,K alpha_i)}  

    The K constants, alpha_1,...,alpha_K, must be positive. The K parameters, 
    theta_1,...,theta_K are nonnegative and sum to 1.
    
    Status:
        Alpha
    """
    __slots__ = 'alpha', '_total', '_mean', 
    
    
    

    def __init__(self, alpha) :
        """
        Args:
            - alpha  -- The parameters of the Dirichlet prior distribution.
                        A vector of non-negative real numbers.  
        """
        # TODO: Check that alphas are positive
        #TODO : what if alpha's not one dimensional?
        self.alpha = asarray(alpha, float64)
        
        self._total = sum(alpha)
        self._mean = None
        
        
    def sample(self) :
        """Return a randomly generated probability vector.
        
        Random samples are generated by sampling K values from gamma
        distributions with parameters a=\alpha_i, b=1, and renormalizing. 
    
        Ref:
            A.M. Law, W.D. Kelton, Simulation Modeling and Analysis (1991).
        Authors:
            Gavin E. Crooks <gec@compbio.berkeley.edu> (2002)
        """
        alpha = self.alpha
        K = len(alpha)
        theta = zeros( (K,), float64)

        for k in range(K):
            theta[k] = random.gammavariate(alpha[k], 1.0) 
        theta /= sum(theta)

        return theta

    def mean(self) :
        if  self._mean ==None:
            self._mean = self.alpha / self._total
        return self._mean

    def covariance(self) :
        alpha = self.alpha
        A = sum(alpha)
        #A2 = A * A
        K = len(alpha)
        cv = zeros( (K,K), float64) 
        
        for i in range(K) :
            cv[i,i] = alpha[i] * (1. - alpha[i]/A) / (A * (A+1.) )
        
        for i in range(K) :
            for j in range(i+1,K) :
                v = - alpha[i] * alpha[j] / (A * A * (A+1.) )
                cv[i,j] = v
                cv[j,i] = v
        return cv
        
    def mean_x(self, x) :
        x = asarray(x, float64)
        if shape(x) != shape(self.alpha) :
            raise ValueError("Argument must be same dimension as Dirichlet")
        return sum( x * self.mean()) 

    def variance_x(self, x) :
        x = asarray(x, float64)
        if shape(x) != shape(self.alpha) :
            raise ValueError("Argument must be same dimension as Dirichlet")        
            
        cv = self.covariance()
        var = na.dot(na.dot(na.transpose( x), cv), x)
        return var

    def odds_ratio(self, composition):
        return [log(value, 2) for value in self.mean() / composition]

    def mean_entropy(self) :
        """Calculate the average entropy of probabilities sampled
        from this Dirichlet distribution. 
        
        Returns:
            The average entropy.
            
        Ref:
            Wolpert & Wolf, PRE 53:6841-6854 (1996) Theorem 7
            (Warning: this paper contains typos.)
        Status:
            Alpha
        Authors:
            GEC 2005
    
        """
        # TODO: Optimize
        alpha = self.alpha
        A = float(sum(alpha))
        ent = 0.0
        for a in alpha:
            if a>0 : ent += - 1.0 * a * digamma( 1.0+a) # FIXME: Check
        ent /= A
        ent += digamma(A+1.0)
        return ent    



    def variance_entropy(self):
        """Calculate the variance of the Dirichlet entropy. 

        Ref:
            Wolpert & Wolf, PRE 53:6841-6854 (1996) Theorem 8
            (Warning: this paper contains typos.)
        """
        alpha = self.alpha
        A = float(sum(alpha))
        A2 = A * (A+1)
        L = len(alpha)
        
        dg1 = zeros( (L) , float64)
        dg2 = zeros( (L) , float64)
        tg2 = zeros( (L) , float64)        
        
        for i in range(L) :
            dg1[i] = digamma(alpha[i] + 1.0)
            dg2[i] = digamma(alpha[i] + 2.0)
            tg2[i] = trigamma(alpha[i] + 2.0)

        dg_Ap2 = digamma( A+2. )
        tg_Ap2 = trigamma( A+2. )
        
        mean = self.mean_entropy()
        var = 0.0
    
        for i in range(L) :
            for j in range(L) :
                if i != j :
                    var += (
                        ( dg1[i] - dg_Ap2 ) * (dg1[j] - dg_Ap2 ) - tg_Ap2 
                        ) * (alpha[i] * alpha[j] ) / A2
                else : 
                    var += (
                        ( dg2[i] - dg_Ap2 ) **2 + ( tg2[i] - tg_Ap2 )
                        ) * ( alpha[i] * (alpha[i]+1.) ) / A2

        var -= mean**2
        return var
        
        
        
    def mean_relative_entropy(self, pvec) :
        ln_p = na.log(pvec)
        return - self.mean_x(ln_p) - self.mean_entropy() 
    
    
    def variance_relative_entropy(self, pvec) :
        ln_p = na.log(pvec)
        return self.variance_x(ln_p) + self.variance_entropy()
    
    
    def interval_relative_entropy(self, pvec, frac) :
        mean = self.mean_relative_entropy(pvec) 
        variance = self.variance_relative_entropy(pvec) 
        sd = sqrt(variance)
        
        # If the variance is small, use the standard 95% 
        # confidence interval: mean +/- 1.96 * sd
        if mean/sd >3.0 :
            return max(0.0, mean - sd*1.96), mean + sd*1.96
        
        g = Gamma.from_mean_variance(mean, variance)
        low_limit = g.inverse_cdf( (1.-frac)/2.)
        high_limit = g.inverse_cdf( 1. - (1.-frac)/2. )
        
        return low_limit, high_limit


class Gamma(object) :
    """The gamma probability distribution. (Not to be confused with the
    gamma function.)
    
    
    """
    __slots__ = 'alpha', 'beta'
    
    def __init__(self, alpha, beta) :
        if alpha <=0.0 :
            raise ValueError("alpha must be positive")
        if beta <=0.0 :
            raise ValueError("beta must be positive")
        self.alpha = alpha
        self.beta = beta
        

    def from_shape_scale(cls, shape, scale) :
        return cls( shape, 1./scale)
    from_shape_scale = classmethod(from_shape_scale)
    

    def from_mean_variance(cls, mean, variance) :
        alpha = mean **2 / variance
        beta = alpha/mean
        return cls(alpha, beta)
    from_mean_variance = classmethod(from_mean_variance)

    def mean(self) :
        return self.alpha / self.beta
        
    def variance(self) :
        return self.alpha / (self.beta**2)

        
    def sample(self) :
        return random.gammavariate(self.alpha, 1./self.beta) 

    def pdf(self, x) :
        if x==0.0 : return 0.0
        a = self.alpha
        b = self.beta
        return (x**(a-1.)) * exp(- b*x )* (b**a) / gamma(a) 
        
    def cdf(self, x) :
        return 1.0-normalized_incomplete_gamma(self.alpha, self.beta*x)

    def inverse_cdf(self, p) :
        def rootof(x) :
            return self.cdf(exp(x)) - p
    
        return exp(find_root(  rootof, log(self.mean()) ) )

#

def find_root(f, x, y=None, fprime=None, tolerance=1.48e-8, max_iterations=50):
    """Return argument 'x' of the function f(x), such that f(x)=0 to
    within the given tolerance. 
    
    f : The function to optimize, f(x)
    x : The initial guess
    y : An optional second guess that shoudl bracket the root.
    fprime : The derivate of f'(x) (Optional)
    tolerance : The error bounds 
    max_iterations : Maximum number of iterations
    
    Raises:
        ArithmeticError :
            Failure to converge to the given tolerence

    Notes:
        Uses Newton-Raphson algorihtm if f'(x) is given, else uses bisect if
        y is given and brackets the root, else uses secant. 

    Status : Beta (Not fully tested)
    """

  
    def secant (f, x, tolerance, max_iterations):
        x0 = x
        x1 = x0+ 1e-4
        v0 = f(x0)
        v1 = f(x1)
        x2 = 0
        
        for i in range(max_iterations):
            # print x0, x1, v0, v1, x2-x0
            x2 = x1 - v1*(x1-x0)/(v1-v0)
            if abs(x2-x1) < tolerance : return x2        
            x0 = x1
            v0 = v1        
            x1 = x2
            v1 = f(x1)

        raise ArithmeticError(
            "Failed to converge after %d iterations, value is %f" \
                % (max_iterations,x1) )


    def bisect(f, a, b, tolerance, max_iterations) :
        fa = f(a)
        fb = f(b)

        if fa==0: return a
        if fb==0: return b
                
        if fa*fb >0 :
            raise ArithmeticError("Start points do not bracket root.")
        
        for i in range(max_iterations):
            #print a,b, fa, fb
            delta = b-a
            xm = a + 0.5*delta # Minimize roundoff in computing the midpoint
            fm = f(xm)
            if delta < tolerance : return xm
            
            if fm * fa >0 : #   Root lies in interval [xm,b], replace a
                a = xm
                fa = fm
            else : #   Root lies in interval [a,xm], replace b
                b = xm
                fb = fm

        raise ArithmeticError(
            "Failed to converge after %d iterations, value is %f" \
                % (max_iterations, x) )
    
    
    def newton(f, x, fprime, tolerance, max_iterations) :
        x0 = x
        for i in range(max_iterations) :
            x1 = x0 - f(x0)/fprime(x0)
            if abs(x1-x0) < tolerance: return x1
            x0 = x1
            
        raise ArithmeticError(
            "Failed to converge after %d iterations, value is %f" \
                % (max_iterations,x1) )

    if fprime is not None :
        return newton(f, x, fprime, tolerance, max_iterations)
    elif y is not None :
        return bisect(f, x, y, tolerance, max_iterations)
    else :
        return secant(f, x, tolerance, max_iterations)
