function rv = rv_coefficient(A, B)
% see  Escoufier, Y. (1973). "Le Traitement des Variables Vectorielles". Biometrics. International Biometric Society. 29 (4): 751–760. doi:10.2307/2529140. JSTOR 2529140. 
% and https://en.wikipedia.org/wiki/RV_coefficient#Definitions
rv = trace(cov(A,B)*cov(B,A))./sqrt(trace(cov(A,A)*cov(A,A))*trace(cov(B,B)*cov(B,B)));