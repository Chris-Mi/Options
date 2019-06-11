import numpy as np
import pandas as pd
from scipy.stats import norm

class vanilla_option:
# This function, developed while learning the basics of object oriented programming, allows to create an option as an object
# with the following inputs:
#    
#  S             1x1   Spot price of the underlying asset of the option at time t = 0
#  K             1x1   Strike price of the option
#  T             1x1   Time to maturity of the option (in years)
#  r             1x1   Risk-free (domestic) interest rate (in decimals, e.g. 0.01 for 1%)
#  q             1x1   Dividend yield / risk-free foreign interest rate (in decimals, e.g. 0.03 for 3%)
#  sd            1x1   Standard deviation (annual) of the underlying asset (in decimals, e.g. 0.15 for 15%)        
#  option_type   Str   Type of the contract: 'Call' or 'Put'
#
# Example:  my_option = vanilla_option(95, 100, 1/12, 0.01, 0.03, 0.18, 'Call')   -> It creates the option object with name "my_option"
#           my_option.K   -> It reminds you about the strike price you used as input
#           my_option.bs_price()   -> It computes the option price using the Black-Scholes formula
    
    
    def __init__(self, S, K, T, r, q, sd, option_type):
        # Initialize object with common option attributes
        self.S = float(S)
        self.K = float(K)
        self.T = float(T)
        self.r = float(r)
        self.q = float(q)
        self.sd = float(sd)
        self.option_type = str(option_type)
        self.d1 = float(( (np.log(self.S/self.K) + (self.r-self.q+(1/2)*self.sd**2)*self.T ) / ( self.sd*np.sqrt(self.T))))
        self.d2 = float(self.d1 - self.sd*np.sqrt(self.T))
        self.moneyness = float(self.S/self.K)
        
        
    def intrinsic(self):
        # Function which computes the intrinsic value of the option
        if self.option_type == 'Call':
            intrinsic = max(self.S-self.K, 0)
        elif self.option_type == 'Put':
            intrinsic = max(self.K-self.S, 0)
        return intrinsic
  
    
    def extrinsic(self):
        # Function which computes the extrinsic (time) value of the option
        extrinsic = self.bs_price() - self.intrinsic()
        return extrinsic
   

    def the_money(self):
        # Function which tells the % moneyness (ITM or OTM) of the option
        if self.option_type == 'Call':
            if self.moneyness == 1:
                the_money = 'The ' + option_type + ' option is At-The-Money (ATM)'
            elif self.moneyness < 1:
                the_money = 'The ' + option_type + ' option is ' + str(round(float((1-self.moneyness)*100),1)) + '% Out-of-The-Money (OTM)'
            elif self.moneyness > 1:
                the_money = 'The ' + option_type + ' option is ' + str(round(float((self.moneyness-1)*100),1)) + '% In-The-Money (ITM)'
        elif self.option_type == 'Put':
            if self.moneyness == 1:
                the_money = 'The ' + option_type + ' option is At-The-Money (ATM)'
            elif self.moneyness < 1:
                the_money = 'The ' + option_type + ' option is ' + str(round(float((1-self.moneyness)*100),1)) + '% In-The-Money (ITM)'
            elif self.moneyness > 1:
                the_money = 'The ' + option_type + ' option is ' + str(round(float((self.moneyness-1)*100),1)) + '% Out-of-The-Money (OTM)'
        return the_money      
        
        
    def bs_price(self):
        # Function which gives the option price according to the Black-Scholes pricing formula
        if self.option_type == 'Call':
            bs_price = self.S*np.exp(-self.q*self.T)*norm.cdf(self.d1) - self.K*np.exp(-self.r*self.T)*norm.cdf(self.d2)
        elif self.option_type == 'Put':
            bs_price = self.K*np.exp(-self.r*self.T)*norm.cdf(-self.d2) - self.S*np.exp(-self.q*self.T)*norm.cdf(-self.d1)
        return bs_price


    def bs_delta(self):
        # Function which gives the delta of the option in the Black-Scholes setting
        if self.option_type == 'Call':
            bs_delta = np.exp(-self.q*self.T)*norm.cdf(self.d1) # delta = N(d1)
        elif self.option_type == 'Put':
            bs_delta = np.exp(-self.q*self.T)*(norm.cdf(self.d1) - 1) # delta(put) = delta(call) - 1
        return bs_delta
       
        
    def bs_gamma(self):
        # Function which gives the gamma of the option in the Black-Scholes setting
        bs_gamma = (np.exp(-self.q*self.T)*norm.pdf(self.d1)) / (self.S*self.sd*np.sqrt(self.T)) # gamma(put) = gamma(call)
        return bs_gamma
     
        
    def bs_vega(self):
        # Function which gives the vega of the option in the Black-Scholes setting
        bs_vega = self.S*np.sqrt(self.T)*np.exp(-self.q*self.T)*norm.pdf(self.d1) # vega(put) = vega(call)
        return bs_vega
    
    
    def bs_rho(self):
        # Function which gives the rho of the option in the Black-Scholes setting
        if self.option_type == 'Call':
            bs_rho = self.T*self.K*np.exp(-self.r*self.T)*norm.cdf(self.d2)
        elif self.option_type == 'Put':
            bs_rho = -self.T*self.K*np.exp(-self.r*self.T)*norm.cdf(-self.d2)
        return bs_rho
     
        
    def bs_theta(self):
        # Function which gives the theta of the option in the Black-Scholes setting
        if self.option_type == 'Call':
            bs_theta = ((-self.S*np.exp(-self.q*self.T)*self.sd*norm.pdf(self.d1) / (2*np.sqrt(self.T)) ) + (self.q*self.S*np.exp(-self.q*self.T)*norm.cdf(self.d1)) / 365 - ( self.r*self.K*np.exp(-self.r*self.T)*norm.cdf(self.d2) )) / 365         
        elif self.option_type == 'Put':
            bs_theta = ((-self.S*self.sd*norm.pdf(self.d1) / (2*np.sqrt(self.T)) ) - (self.q*self.S*np.exp(-self.q*self.T)*norm.cdf(self.d1)) / 365 + ( self.r*self.K*np.exp(-self.r*self.T)*norm.cdf(-self.d2) )) / 365
        return bs_theta
       
        
    def mc_price(self, number_sim, method):
        # Function which gives the price of the option using Monte Carlo simulation and the additional inputs:
#         number_sim   1x1   Number of simulations
#         method       Str   Method used for efficient simulation, either 'Naive' (simple simulation) or 'Anti' (more efficient
#                            simulation method using antithetic variates)
#
#         Example:  my_option = vanilla_option(95, 100, 1/12, 0.01, 0.03, 0.18, 'Call')
#                   my_option.mc_price(10000, 'Anti')
#                 or directly:
#                   vanilla_option(95, 100, 1/12, 0.01, 0.03, 0.18, 'Call').mc_price(10000, 'Anti')
        drift = (self.r - self.q - 0.5*self.sd**2)*self.T
        vola = self.sd * np.sqrt(self.T)
        epsilons = pd.Series(np.random.normal(loc = 0, scale = 1, size = number_sim))      
        if method == 'Anti':
            sim_prices = pd.Series([self.S*np.exp(drift + vola*eps) for eps in epsilons])
            sim_prices_anti = pd.Series([self.S*np.exp(drift + vola*(-eps)) for eps in epsilons])
            if self.option_type == 'Call':
                payoffs = pd.Series([max(payoff, 0) for payoff in (sim_prices - self.K)])
                payoffs_anti = pd.Series([max(payoff, 0) for payoff in (sim_prices_anti - self.K)])
                mc_price = (0.5*payoffs + 0.5*payoffs_anti).mean()
            elif self.option_type == 'Put':
                payoffs = pd.Series([max(payoff, 0) for payoff in (self.K - sim_prices)])
                payoffs_anti = pd.Series([max(payoff, 0) for payoff in (self.K - sim_prices_anti)])
                mc_price = (0.5*payoffs + 0.5*payoffs_anti).mean()                    
        elif method == 'Naive':
            sim_prices = pd.Series([self.S*np.exp(drift + vola*eps) for eps in epsilons])
            if self.option_type == 'Call':
                mc_price = pd.Series([max(payoff, 0) for payoff in (sim_prices - self.K)]).mean()
            elif self.option_type == 'Put':
                mc_price = pd.Series([max(payoff, 0) for payoff in (self.K - sim_prices)]).mean()        
        return mc_price    