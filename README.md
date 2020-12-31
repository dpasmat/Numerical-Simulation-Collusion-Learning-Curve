# Numerical-Simulation-Collusion-Learning-Curve

Learning-by-doing is the process whereby a firm's *current* output reduces its *future* cost. Learning-by-doing is an important feature of the production process in many computer and electronics products. Moreover, in recent years, numerous manufacturers of these products have been subject to fines for explicitly colluding to raise prices.  

This README describes how to run the numerical simulation model in the paper “Collusion along the learning curve: Theory and evidence from the semiconductor industry.” The model calculates the minimum discount factor necessary to sustain collusion in a duopoly model, as a function of firms' rate of learning-by-doing. 

This output from this code produces the data underlying Figure 8 from the paper. A modified version of this code can also be used to produce Figure 9(a). However, this modification is not described in this README and is available upon request. 

# Description of programs

The program `Solve_model.m` specifies and solves the linear demand version of the repeated game model of collusion that is constructed in the paper. The program contains the supply and demand parameters of the model, to be specified by the user.

- There are two loops in the script: the outer loop is over the variable `dt` – the discount factor. The inner loop is over the variable `ll` - which specifies whether the firm is adhering to the collusive path or deviating; and if deviating, for which product the firm is deviating.

- To solve the model, the script specifies the first order conditions for each firm, product, and time period. It then uses the function `Solve` to find the vector of equilibrium outputs. Lastly, for each iteration of the outer loop `dt`, it computes the difference in profitability between staying on the equilibrium path and deviating, separately for each product. These values are stored as elements in each row of the matrix `comp`.

The program `profit_function.m` calculates the market prices for each product and time period; and the firm’s profits for each time period. It must be saved in the same directory as `Solve_model.m`, as in this repository.

# Software Requirements

Matlab (code was last run with Version R2020a)
- Symbolic Math Toolbox
	
# Memory and Runtime Requirements

The code was last run on a dual processor, 4-core Intel Xeon server with 64Gb of RAM at 3.0 GHz, running on Red Hat Linux 5.3. Computation took 4.25 hours.

# Other Notes

The first order conditions (FOCs) in `Solve_model.m` can be changed to accomodate different assumptions on the demand or supply side of the economic model.
- On the supply side, the FOCs can be updated to reflect a multiproduct (rather than single product) punishment strategy. The equations are represented in Web Appendix B of the paper. 
- On the demand side, the FOCs can be updated to reflect different functional forms. To guarantee the existence of an equilibrium, a sufficient condition is log-concave inverse demand, e.g. with the exponential or log-linear form. Further background is outlined in Web Appendix C of the paper. 



