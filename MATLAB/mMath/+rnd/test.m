classdef test < matlab.unittest.TestCase
  %test Unit tests for methods in rnd package
  %   Tests include:
  %   - Monte Carlo test of validity of samples
  
  properties
  end
  
  methods (Test)
    
    function test_gauss( test )
      
      % Take a random mean and covariance
      dim = 6;
      muVec = randn(dim,1);
      covMat = rnd.cov(dim);
%       covMat = eye(3);
      
      numSamples = 1e6;
      V = rnd.gauss(muVec,covMat,numSamples);
      % Compute sample mean and covariance
      m = 1/numSamples* sum(V,2);
      V_m = V-repmat(m,1,numSamples);
      c = 1/(numSamples-1)* (V_m*V_m');
%       m = mean(V')';
%       c = cov(V');
      
      % Compare the parameters from population to the original parameters
%       m - muVec
%       c - covMat
%       keyboard

      test.verifyEqual(norm(m-muVec),0,'AbsTol',1e-2);
      test.verifyEqual(norm(c-covMat),0,'AbsTol',1e-2);
      
    end
  
  end
  
end