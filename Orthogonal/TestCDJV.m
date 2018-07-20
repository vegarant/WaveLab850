classdef TestCDJV < matlab.unittest.TestCase
    % TestCDJV test if a the DWT and IDWT function calles reveres each other
    % for all kernel functions. It also checks if the vanishing moments 
    % condition holds for preconditioned wavelets.
    %
    % USAGE:
    % >> testCase = TestCDJV;
    % >> res = run(testCase)
    %
    properties
        x, N, eps, nres, dim, J0;
    end  
    
    methods (Test)

        % Constructor
        function obj = TestCDJV(testCase)
            obj.dim = 7;
            obj.N = 2^obj.dim;
            obj.x = rand([obj.N,1]);
            obj.eps = 1e-5;
            obj.nres = 2;
            obj.J0 = obj.dim - obj.nres;
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function testDB2_bd(testCase)
            testWaveletDecRec(testCase, 2);
        end
        function testDB3_bd(testCase)
            testWaveletDecRec(testCase, 3);
        end
        function testDB4_bd(testCase)
            testWaveletDecRec(testCase, 4);
        end
        function testDB5_bd(testCase)
            testWaveletDecRec(testCase, 5);
        end
        function testDB6_bd(testCase)
            testWaveletDecRec(testCase, 6);
        end
        function testDB7_bd(testCase)
            testWaveletDecRec(testCase, 7);
        end
        function testDB8_bd(testCase)
            testWaveletDecRec(testCase, 8);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function testDB2_bd_noP(testCase)
            testWaveletDecRec_noP(testCase, 2);
        end
        function testDB3_bd_noP(testCase)
            testWaveletDecRec_noP(testCase, 3);
        end
        function testDB4_bd_noP(testCase)
            testWaveletDecRec_noP(testCase, 4);
        end
        function testDB5_bd_noP(testCase)
            testWaveletDecRec_noP(testCase, 5);
        end
        function testDB6_bd_noP(testCase)
            testWaveletDecRec_noP(testCase, 6);
        end
        function testDB7_bd_noP(testCase)
            testWaveletDecRec_noP(testCase, 7);
        end
        function testDB8_bd_noP(testCase)
            testWaveletDecRec_noP(testCase, 8);
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function testDB2_boundary_property(testCase)
            testPrecondBoundary(testCase, 2);  
        end
        function testDB3_boundary_property(testCase)
            testPrecondBoundary(testCase, 3);  
        end
        function testDB4_boundary_property(testCase)
            testPrecondBoundary(testCase, 4);  
        end
        function testDB5_boundary_property(testCase)
            testPrecondBoundary(testCase, 5);  
        end
        function testDB6_boundary_property(testCase)
            testPrecondBoundary(testCase, 6);  
        end
        function testDB7_boundary_property(testCase)
            testPrecondBoundary(testCase, 7);  
        end
        function testDB8_boundary_property(testCase)
            testPrecondBoundary(testCase, 8);  
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function testDB2_orthonormality(testCase)
            testOrthonormality(testCase, 2);
        end
        function testDB3_orthonormality(testCase)
            testOrthonormality(testCase, 3);
        end
        function testDB4_orthonormality(testCase)
            testOrthonormality(testCase, 4);
        end
        function testDB5_orthonormality(testCase)
            testOrthonormality(testCase, 5);
        end
        function testDB6_orthonormality(testCase)
            testOrthonormality(testCase, 6);
        end
        function testDB7_orthonormality(testCase)
            testOrthonormality(testCase, 7);
        end
        function testDB8_orthonormality(testCase)
            testOrthonormality(testCase, 8);
        end

    end

    methods (Access=private)
        
        
        function testWaveletDecRec(testCase, vm)
            x   = testCase.x;
            eps = testCase.eps;
            J0  = testCase.J0;
            z = IWT_CDJV( FWT_CDJV(x, J0, vm), J0, vm);
            % Check relative error
            testCase.verifyTrue(norm(z-x,2) < norm(x,2)*eps);
        end

        function testWaveletDecRec_noP(testCase, vm)
            x   = testCase.x;
            eps = testCase.eps;
            J0  = testCase.J0;
            z = IWT_CDJV_noP( FWT_CDJV_noP(x, J0, vm), J0, vm);
            % Check relative error
            testCase.verifyTrue(norm(z-x,2) < norm(x,2)*eps);
        end

        function testPrecondBoundary(testCase, vm)
            
            eps = testCase.eps;
            nres = testCase.nres;
            J0 = testCase.J0;
            N = testCase.N;

            x = linspace(0,1,N)';
            success = 1;
            for i = 1:vm
                y = x.^(i-1);
                z = FWT_CDJV(x, J0, vm);
                a = z(N/2^nres+1:N);
                %fprintf('vm: %d, i: %d, ||a||: %g\n', vm, i, norm(a,2));
                success = success & norm(a,2) < eps;
            end
            testCase.verifyTrue(success);
            
        end
        
        function testOrthonormality(testCase, vm)
            N   = testCase.N;
            J0  = testCase.J0;
            eps = testCase.eps; 
            
            W = zeros([N,N]);
            for i = 1:N
                ei = zeros([N,1]); ei(i)=1;
                W(:,i) = FWT_CDJV_noP(ei, J0, vm);
            end
            A = W'*W;
            B = W*W';
            I = eye(N);
            %disp('-----------')
            %vm
            %norm(A-I, 'fro')
            %norm(B-I, 'fro')
            testCase.verifyTrue(norm(A-I,'fro') < eps);
            testCase.verifyTrue(norm(B-I,'fro') < eps);
        end

    end
end 



