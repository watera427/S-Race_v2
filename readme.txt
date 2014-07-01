S_Race_v2 is a simple implementation of a multi-objective racing procesure based on sign test, discrete holm's step-down procedure, and an alpha adapting scheme.

For detail, please refer to the following paper:

T. Zhang, M. Georgiopoulos, G. C. Anagnostopoulos, "Multi-Objective Model Selection via Racing", 2014 

S_Race_v2 is available at https://github.com/watera427/S-Race_v2

Author contact: Tiantian Zhang 

Email: zhangtt@knights.ucf.edu

COPYRIGHT (C) 2014 Tiantian Zhang, Machine Learning Lab, University of Central Florida
All rights reserved.

Redistribution and use in source and binary forms, with or without 
modification, are permitted provided that the following conditions are 
met:

+ Redistributions of source code must retain the above copyright 
  notice, this list of conditions and the following disclaimer.
+ Redistributions in binary form must reproduce the above copyright 
  notice, this list of conditions and the following disclaimer in 
  the documentation and/or other materials provided with the distribution
      
THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE 
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF 
SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
POSSIBILITY OF SUCH DAMAGE.

=======================================================================================



S_Race_v2 methods:
================

+ obj = S_Race_v2(M, Max_step, Delta, batch_size)

  Create a S_Race_v2 object before starting a racing

  ' M ' is the number of initial models
  ' Max_step ' is the predefined maximal number of steps in racing
  ' Delta ' is the predefined accuracy of S-Race
  ' batch_size' is the number of test instances for each step

+ retained = Racing(obj, results)
   
  Every step, call 'Racing' with newly generated results. The index of retained models will   
  be returned 

  'results' should be stored as a 3-dimensional matrix of size batch_size x M x no_obj
  'batch_size' is the number of test instance determined by user
  'M' is the number of retained models from last step. The initial M equals to the number of initial models
  'no_obj' is the number of objectives considered when comparing models



S_Race variables:
==================
        M               % number of current models
        currentM        % index of current models, initialized as 1:M        
        batch_size      % number of the test instances every step
        results         % a 3 dimemsional matrix of size (step x batch_size) x M x No_obj matrix, storing the performance vectors of all current models so far  
        Max_step        % assign the maximum number of steps as the stopping criterion
        alpha           % alpha value used in each family of S-Race
        No_obj          % number of objectives
        Delta           % overall probability of making any Type I errors
        fam             % the number of families of current step



Tips on Practical Use:
==================
+ Please organize the results in a 3-D matrix as required
+ S_Race_v2 is designed for minimization problem only
+ Do not reinitalize S_Race_v2 object before racing is complete
+ Remember to add path in default search path ('STARTUP.m')
+ Author used MATLAB R2013b. 



Examples:
==================
Please refer to 'example.m'





